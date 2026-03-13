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
#include <deque>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>


// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only DNA chunk (no N, no newlines — pre-filtered by helicase),
// detect superkmer boundaries with the minimizer iterator, and append each
// superkmer to the corresponding writer.
//
// Splits on every minimizer HASH change (not merely partition change) so that
// every superkmer has a single well-defined minimizer whose position can be
// stored in the 1-byte min_pos header field.  Superkmers that cross a minimizer
// boundary but stay in the same partition are thus split into two shorter pieces;
// each piece is routed to the same partition file — correctness is unchanged.
//
// find_minimizer_pos: O(superkmer_len) scan returning the 0-indexed start of
// the leftmost l-mer with minimum canonical ntHash in seq[0..len-1].

template <uint16_t l>
static uint8_t find_minimizer_pos(const char* seq, uint8_t len) noexcept
{
    nt_hash::Roller<l> roller;
    roller.init(seq);
    uint64_t min_h = roller.canonical();
    uint8_t  min_p = 0;
    for (uint8_t i = l; i < len; ++i) {
        roller.roll(nt_hash::to_2bit(seq[i - l]), nt_hash::to_2bit(seq[i]));
        const uint64_t h = roller.canonical();
        if (h < min_h) { min_h = h; min_p = static_cast<uint8_t>(i - l + 1); }
    }
    return min_p;
}

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
    uint64_t prev_hash = min_it.hash();
    size_t   pid       = partition_fn(prev_hash);
    size_t   sk_start  = 0;
    size_t   sk_end    = k;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance(seq[pos]);
        const uint64_t new_hash = min_it.hash();
        // Flush on minimizer change OR when sk_len is about to exceed uint8_t max.
        // At flush time sk_end == pos, so sk_len = pos - sk_start ≤ 255 always.
        if (new_hash != prev_hash || pos - sk_start >= 255u) {
            const auto sk_len = static_cast<uint8_t>(sk_end - sk_start);
            writers[pid].append(seq + sk_start, sk_len,
                                find_minimizer_pos<l>(seq + sk_start, sk_len));
            kmer_count += sk_len - k + 1;
            prev_hash = new_hash;
            pid       = partition_fn(new_hash);
            sk_start  = pos - (k - 1);
        }
        sk_end = pos + 1;
    }

    const auto sk_len = static_cast<uint8_t>(sk_end - sk_start);
    writers[pid].append(seq + sk_start, sk_len,
                        find_minimizer_pos<l>(seq + sk_start, sk_len));
    kmer_count += sk_len - k + 1;
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// Producer/consumer design: the calling thread reads all input files and pushes
// ACTG-only chunks onto a bounded queue; num_threads workers drain the queue in
// parallel.  Each worker owns its MinimizerWindow and SuperkmerWriters, flushing
// to the shared bucket ofstreams under per-bucket mutexes.
//
// Chunks are copied into the queue so their lifetime is independent of the
// parser / mmap.  For a 3 GB human genome the copy adds ~300 ms at ~10 GB/s
// memcpy bandwidth — negligible against the processing time.  The queue cap
// (QUEUE_CAP tasks) bounds peak memory to a few hundred MB at most.

template <uint16_t k, uint16_t l, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn)
{
    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts   = cfg.num_partitions;

    // ── Work queue ────────────────────────────────────────────────────
    static constexpr size_t QUEUE_CAP = 128;
    struct Task { std::vector<char> data; };

    std::mutex              q_mtx;
    std::condition_variable q_not_empty, q_not_full;
    std::deque<std::unique_ptr<Task>> q;
    bool q_done = false;

    std::vector<std::mutex> bucket_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    // ── Workers ───────────────────────────────────────────────────────
    auto worker = [&]() {
        MinimizerWindow<k, l>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts);
        uint64_t local_seqs = 0, local_kmers = 0;

        while (true) {
            std::unique_ptr<Task> task;
            {
                std::unique_lock<std::mutex> lk(q_mtx);
                q_not_empty.wait(lk, [&]{ return !q.empty() || q_done; });
                if (q.empty()) break;
                task = std::move(q.front());
                q.pop_front();
            }
            q_not_full.notify_one();

            extract_superkmers_from_actg<k, l>(
                task->data.data(), task->data.size(),
                partition_fn, min_it, writers, local_kmers);
            ++local_seqs;

            for (size_t p = 0; p < n_parts; ++p)
                if (writers[p].needs_flush())
                    writers[p].flush_to(buckets[p], bucket_mutexes[p]);
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

    // ── Producer (calling thread) ──────────────────────────────────────
    SeqSource source;
    for (const auto& path : cfg.input_files) {
        source.process(path, [&](const char* chunk, size_t len) {
            auto task = std::make_unique<Task>();
            task->data.assign(chunk, chunk + len);
            {
                std::unique_lock<std::mutex> lk(q_mtx);
                q_not_full.wait(lk, [&]{ return q.size() < QUEUE_CAP; });
                q.push_back(std::move(task));
            }
            q_not_empty.notify_one();
        });
    }

    {
        std::lock_guard<std::mutex> lk(q_mtx);
        q_done = true;
    }
    q_not_empty.notify_all();

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
