#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// This file implements the two strategies that share the same minimizer-hash
// inner loop, differing only in how they map a minimizer hash to a partition ID:
//
//   partition_kmers_hash<k,l>(cfg, buckets)           — hash % num_partitions
//   partition_kmers_kmtricks<k,l>(cfg, buckets, rt)   — RepartitionTable lookup
//
// Parsing is done with helicase (SIMD FASTA/FASTQ, ~5-7 GB/s).
// helicase with SPLIT_NON_ACTG delivers ACTG-only chunks (N-free, newline-free)
// straight to extract_superkmers_from_actg — no placeholder checks needed.
//
// Shared internals:
//   extract_superkmers_from_actg<k, PartitionFn>  — pure ACTG sequence logic,
//                                                    no I/O, no threads.
//   partition_kmers_impl<k, l, PartitionFn>        — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "partition_kmtricks.hpp"
#include "minimizer_window.hpp"

#include <helicase.hpp>
#include "fast_fasta.hpp"

#include <zlib.h>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <cstring>


struct PartitionStats { uint64_t seqs = 0, kmers = 0; };


// ─── helicase config ──────────────────────────────────────────────────────────
//
// DNA string + split at non-ACTG + return each chunk separately.
// One next() call per ACTG-only run (newline-free, N-free).

static constexpr helicase::Config HELICASE_ACTG =
    helicase::ParserOptions()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .config()
    & ~helicase::advanced::RETURN_RECORD;


// FileBuffer removed: gz files are now streamed via SeqReader (64 MB window)
// to avoid OOM on large compressed inputs (e.g. 5.9 GB gz → ~20 GB decompressed).


// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only DNA chunk (no N, no newlines — pre-filtered by helicase),
// detect superkmer boundaries with the minimizer iterator, and append each
// superkmer to the corresponding writer.
//
// Because the input is guaranteed ACTG, there is no placeholder check and no
// need for the in_run / skip-k-ahead startup; the window opens immediately.

template <uint16_t k, typename PartitionFn>
void extract_superkmers_from_actg(
    const char* const             seq,
    const size_t                  seq_len,
    PartitionFn&&                 partition_fn,
    MinimizerWindow<k>&           min_it,
    std::vector<SuperkmerWriter>& writers,
    uint64_t&                     kmer_count)
{
    if (seq_len < k) return;

    min_it.reset(seq);
    size_t pid     = partition_fn(min_it.hash());
    size_t sk_start = 0;
    size_t sk_end   = k;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance(seq[pos]);
        const size_t new_pid = partition_fn(min_it.hash());
        if (new_pid != pid) {
            writers[pid].append(seq + sk_start,
                                static_cast<uint32_t>(sk_end - sk_start));
            kmer_count += sk_end - sk_start - k + 1;
            pid      = new_pid;
            sk_start = pos - (k - 1);
        }
        sk_end = pos + 1;
    }

    writers[pid].append(seq + sk_start, static_cast<uint32_t>(sk_end - sk_start));
    kmer_count += sk_end - sk_start - k + 1;
}


// ─── Legacy: placeholder-aware version (kept for reference / non-helicase paths) ─

template <uint16_t k, typename PartitionFn>
void extract_superkmers_from_seq(
    const char* const             seq,
    const size_t                  seq_len,
    PartitionFn&&                 partition_fn,
    MinimizerWindow<k>&           min_it,
    std::vector<SuperkmerWriter>& writers,
    uint64_t&                     kmer_count)
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

        if (!nt_hash::is_dna(ch)) {
            if (in_run) { emit_superkmer(); in_run = false; }
            continue;
        }

        if (!in_run) {
            if (pos + k > seq_len) break;

            bool ok = true;
            for (size_t t = 1; t < k; ++t) {
                if (!nt_hash::is_dna(seq[pos + t])) {
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


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// Spawns min(num_threads, n_files) threads. Each thread processes its assigned
// input files round-robin. helicase delivers ACTG-only chunks via SliceInput;
// extract_superkmers_from_actg maps each chunk to superkmers + partition writers.

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
        MinimizerWindow<k> min_it(l);
        std::vector<SuperkmerWriter> writers(cfg.num_partitions);
        SeqReader sreader; // reused across gz files; only grows, never shrinks

        uint64_t local_seqs = 0, local_kmers = 0;

        const auto flush_if_needed = [&]() {
            for (size_t p = 0; p < cfg.num_partitions; ++p)
                if (writers[p].needs_flush())
                    writers[p].flush_to(buckets[p], bucket_mutexes[p]);
        };

        for (size_t fi = tid; fi < n_files; fi += n_threads) {
            const auto& path = cfg.input_files[fi];
            const bool is_gz = path.size() > 3 &&
                               path.compare(path.size() - 3, 3, ".gz") == 0;

            if (is_gz) {
                // Streaming gz: SeqReader (64 MB sliding window, O(1) memory).
                // Use placeholder-aware extract_superkmers_from_seq to handle N's.
                if (!sreader.load(path)) continue;
                while (sreader.read_next_seq()) {
                    ++local_seqs;
                    extract_superkmers_from_seq<k>(
                        sreader.seq(), sreader.seq_len(),
                        partition_fn, min_it, writers, local_kmers);
                    flush_if_needed();
                }
            } else {
                // Plain file: zero-copy mmap + helicase SIMD (ACTG-only chunks).
                helicase::MmapInput inp(path);
                const uint8_t fb = inp.first_byte();
                if (fb == '@') {
                    helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    while (p.next()) {
                        auto [ptr, len] = p.get_dna_raw();
                        ++local_seqs;
                        extract_superkmers_from_actg<k>(
                            ptr, len, partition_fn, min_it, writers, local_kmers);
                        flush_if_needed();
                    }
                } else {
                    helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    while (p.next()) {
                        auto [ptr, len] = p.get_dna_raw();
                        ++local_seqs;
                        extract_superkmers_from_actg<k>(
                            ptr, len, partition_fn, min_it, writers, local_kmers);
                        flush_if_needed();
                    }
                }
            }
        }

        // Final flush of any remaining buffered data.
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
