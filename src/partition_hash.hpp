#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// This file implements the two strategies that share the same minimizer-hash
// inner loop (extract_superkmers_from_seq), differing only in how they map a
// minimizer hash to a partition ID:
//
//   partition_kmers<k, l>(cfg, buckets)          — hash strategy:
//                                                   hash % num_partitions
//   partition_kmers_kmtricks<k, l>(cfg, bkts, rt)— kmtricks strategy:
//                                                   RepartitionTable lookup
//
// Shared internals:
//   extract_superkmers_from_seq<k, PartitionFn>  — pure sequence logic, no I/O,
//                                                   no threads. Takes any callable
//                                                   (uint64_t hash → size_t pid).
//   partition_kmers_impl<k, l, PartitionFn>      — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "partition_kmtricks.hpp"
#include "fast_fasta.hpp"
#include "DNA_Utility.hpp"
#include "Minimizer_Iterator.hpp"

#include <vector>
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>


struct PartitionStats { uint64_t seqs = 0, kmers = 0; };


// ─── Partition logic brick ────────────────────────────────────────────────────
//
// Walk one DNA sequence, detect superkmer boundaries using the minimizer
// iterator, and append each superkmer to the corresponding per-bucket writer.
//
// PartitionFn: any callable with signature  size_t fn(uint64_t hash).
//   Hash strategy:     [n](uint64_t h){ return h % n; }
//   Kmtricks strategy: [&rt](uint64_t h){ return rt(h); }
//
// No I/O, no locking, no thread state: fully self-contained and testable.

template <uint16_t k, typename PartitionFn>
void extract_superkmers_from_seq(
    const char* const               seq,
    const size_t                    seq_len,
    PartitionFn&&                   partition_fn,
    cuttlefish::Min_Iterator<k>&    min_it,
    std::vector<SuperkmerWriter>&   writers,
    uint64_t&                       kmer_count)
{
    if (seq_len < k) return;

    bool   in_run   = false;
    size_t sk_start = 0, sk_end = 0;
    size_t pid      = 0;

    const auto emit_superkmer = [&]() {
        const auto len = static_cast<uint32_t>(sk_end - sk_start);
        writers[pid].append(seq + sk_start, len);
        kmer_count += len - k + 1;
    };

    for (size_t pos = 0; pos < seq_len; ++pos) {
        const char ch = seq[pos];

        if (cuttlefish::DNA_Utility::is_placeholder(ch)) {
            if (in_run) { emit_superkmer(); in_run = false; }
            continue;
        }

        if (!in_run) {
            if (pos + k > seq_len) break;

            bool ok = true;
            for (size_t t = 1; t < k; ++t) {
                if (cuttlefish::DNA_Utility::is_placeholder(seq[pos + t])) {
                    pos += t; ok = false; break;
                }
            }
            if (!ok) continue;

            min_it.reset(seq + pos);
            pid      = partition_fn(min_it.hash());
            sk_start = pos;
            sk_end   = pos + k;
            in_run   = true;
            pos     += k - 1;
            continue;
        }

        min_it.advance(ch);
        const size_t new_pid = partition_fn(min_it.hash());

        if (new_pid != pid) {
            emit_superkmer();
            pid      = new_pid;
            sk_start = pos - (k - 1);
        }
        sk_end = pos + 1;
    }

    if (in_run) emit_superkmer();
}


// ─── Parallel harness (shared implementation) ─────────────────────────────────
//
// Spawns min(num_threads, n_files) threads.  Each thread processes its assigned
// input files round-robin, calling extract_superkmers_from_seq per sequence.
// Per-bucket SuperkmerWriters buffer data locally and flush to the shared
// ofstreams under per-bucket mutexes.
//
// partition_fn: any callable uint64_t → size_t.

template <uint16_t k, uint16_t l, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn)
{
    const size_t n_files   = cfg.input_files.size();
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_files);

    std::vector<std::mutex> bucket_mutexes(cfg.num_partitions);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    auto worker = [&](size_t tid) {
        cuttlefish::Min_Iterator<k> min_it(l);
        std::vector<SuperkmerWriter> writers(cfg.num_partitions);

        uint64_t local_seqs = 0, local_kmers = 0;

        SeqReader parser;
        for (size_t fi = tid; fi < n_files; fi += n_threads) {
            if (!parser.load(cfg.input_files[fi])) continue;

            while (parser.read_next_seq()) {
                ++local_seqs;

                extract_superkmers_from_seq<k>(
                    parser.seq(), parser.seq_len(),
                    partition_fn, min_it, writers, local_kmers);

                for (size_t p = 0; p < cfg.num_partitions; ++p)
                    if (writers[p].needs_flush())
                        writers[p].flush_to(buckets[p], bucket_mutexes[p]);
            }
        }

        for (size_t p = 0; p < cfg.num_partitions; ++p)
            writers[p].flush_to(buckets[p], bucket_mutexes[p]);

        total_seqs .fetch_add(local_seqs,  std::memory_order_relaxed);
        total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    return { total_seqs.load(), total_kmers.load() };
}


// ─── Public: hash strategy ────────────────────────────────────────────────────
//
// Partition by hash % num_partitions — no pre-scan required.

template <uint16_t k, uint16_t l>
PartitionStats partition_kmers_hash(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets)
{
    const size_t n = cfg.num_partitions;
    return partition_kmers_impl<k, l>(cfg, buckets,
        [n](uint64_t h) -> size_t { return h % n; });
}


// ─── Public: kmtricks strategy ────────────────────────────────────────────────
//
// Partition by repartition table lookup — requires a pre-built RepartitionTable
// (either RepartitionTable::uniform() or built from scan_minimizer_counts()).

template <uint16_t k, uint16_t l>
PartitionStats partition_kmers_kmtricks(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    const RepartitionTable&     rt)
{
    return partition_kmers_impl<k, l>(cfg, buckets,
        [&rt](uint64_t h) -> size_t { return rt(h); });
}
