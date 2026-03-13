#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// This file implements the two strategies that share the same minimizer-hash
// inner loop, differing only in how they map a minimizer hash to a partition ID:
//
//   partition_kmers_hash<k,l>(cfg, buckets)           — hash % num_partitions
//   partition_kmers_kmtricks<k,l>(cfg, buckets, rt)   — RepartitionTable lookup
//
// Parsing is delegated to SeqSource, which delivers ACTG-only chunks regardless
// of file format or compression (helicase for plain files, SeqReader+split for gz).
//
// Shared internals:
//   extract_superkmers_from_actg<k, PartitionFn>  — pure ACTG sequence logic,
//                                                    no I/O, no threads.
//   partition_kmers_impl<k, l, PartitionFn>        — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "partition_kmtricks.hpp"
#include "minimizer_window.hpp"
#include "seq_source.hpp"

#include <vector>
#include <thread>
#include <mutex>
#include <atomic>


// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only DNA chunk (no N, no newlines — pre-filtered by helicase),
// detect superkmer boundaries with the minimizer iterator, and append each
// superkmer to the corresponding writer.
//
// Splits on every minimizer HASH change (not merely partition change) so that
// every superkmer has a single well-defined minimizer.
//
// min_pos is obtained from MinimizerWindow::min_lmer_pos() — the position of the
// minimizer l-mer is tracked as a side-effect of the sliding-window minimum with
// no extra scan.  This replaces the old O(sk_len) find_minimizer_pos rescan that
// consumed ~22% of phase1 cycles.

template <uint16_t k, uint16_t l, typename PartitionFn>
void extract_superkmers_from_actg(
    const char* const             seq,
    const size_t                  seq_len,
    PartitionFn&&                 partition_fn,
    MinimizerWindow<k, l>&        min_it,
    std::vector<SuperkmerWriter>& writers,
    uint64_t&                     kmer_count)
{
    if (seq_len < k) return;

    min_it.reset(seq);
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = min_it.min_lmer_pos(); // absolute pos within seq
    size_t   pid          = partition_fn(prev_hash);
    size_t   sk_start     = 0;
    size_t   sk_end       = k;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance(seq[pos]);
        const uint64_t new_hash = min_it.hash();
        // Flush on minimizer change OR when sk_len is about to exceed uint8_t max.
        // At flush time sk_end == pos, so sk_len = pos - sk_start ≤ 255 always.
        if (new_hash != prev_hash || pos - sk_start >= 255u) {
            const auto sk_len  = static_cast<uint8_t>(sk_end - sk_start);
            const auto min_pos = static_cast<uint8_t>(prev_min_pos - sk_start);
            writers[pid].append(seq + sk_start, sk_len, min_pos);
            kmer_count += sk_len - k + 1;
            prev_hash    = new_hash;
            prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
        sk_end = pos + 1;
    }

    const auto sk_len  = static_cast<uint8_t>(sk_end - sk_start);
    const auto min_pos = static_cast<uint8_t>(prev_min_pos - sk_start);
    writers[pid].append(seq + sk_start, sk_len, min_pos);
    kmer_count += sk_len - k + 1;
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// File-level work stealing: each worker atomically claims the next unprocessed
// file and runs helicase + superkmer extraction entirely on its own thread.
// No producer thread, no queue, no memcpy — exactly n_threads threads run,
// all land on P-cores (avoids the hybrid-CPU E-core penalty that the old
// producer/consumer design suffered from when n_threads+1 threads competed
// for n_threads P-cores).
//
// Load balancing is coarse (file granularity) but for equal-size files
// (e.g. 200 × 4.7 MB E. coli) the imbalance is negligible.

template <uint16_t k, uint16_t l, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn)
{
    const size_t n_files   = cfg.input_files.size();
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_files);
    const size_t n_parts   = cfg.num_partitions;

    std::atomic<size_t>   next_file{0};
    std::vector<std::mutex> bucket_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    auto worker = [&]() {
        SeqSource            source;
        MinimizerWindow<k, l>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts);
        uint64_t local_seqs = 0, local_kmers = 0;

        while (true) {
            const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
            if (fi >= n_files) break;

            source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                extract_superkmers_from_actg<k, l>(
                    chunk, len, partition_fn, min_it, writers, local_kmers);
                ++local_seqs;
                for (size_t p = 0; p < n_parts; ++p)
                    if (writers[p].needs_flush())
                        writers[p].flush_to(buckets[p], bucket_mutexes[p]);
            });
        }

        for (size_t p = 0; p < n_parts; ++p)
            writers[p].flush_to(buckets[p], bucket_mutexes[p]);

        total_seqs .fetch_add(local_seqs,  std::memory_order_relaxed);
        total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    return { total_seqs.load(), total_kmers.load() };
}


// ─── Public: hash strategy ────────────────────────────────────────────────────

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

template <uint16_t k, uint16_t l>
PartitionStats partition_kmers_kmtricks(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    const RepartitionTable&     rt)
{
    return partition_kmers_impl<k, l>(cfg, buckets,
        [&rt](uint64_t h) -> size_t { return rt(h); });
}
