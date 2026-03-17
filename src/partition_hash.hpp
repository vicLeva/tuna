#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// Partitioning: minimizer_hash % num_partitions.
//
// Parsing is delegated to SeqSource, which delivers ACTG-only chunks regardless
// of file format or compression (helicase for plain files, SeqReader+split for gz).
//
// Internals:
//   extract_superkmers_from_actg<k, l, PartitionFn>  — pure ACTG sequence logic,
//                                                       no I/O, no threads.
//   partition_kmers_impl<k, l, PartitionFn>           — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "minimizer_window.hpp"
#include "seq_source.hpp"

#include <vector>
#include <deque>
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


// ─── Producer-consumer harness for a single gz file ──────────────────────────
//
// When there is only one input file and it is gzip-compressed, the standard
// work-stealing harness wastes all but one thread (min(n_threads,1)=1).
// This harness instead:
//   • 1 producer thread: gz decompression + split_actg → batches of ACTG chunks
//   • n_threads-1 consumer threads: extract_superkmers_from_actg per chunk
//
// The queue is bounded (MAX_QUEUE batches) for backpressure.  Consumers flush
// SuperkmerWriters to the shared bucket files under per-bucket mutexes.

template <uint16_t k, uint16_t l, typename PartitionFn>
PartitionStats partition_kmers_gz_pc(
    const Config&               cfg,
    const std::string&          gz_path,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      n_threads)      // ≥ 2 (1 producer + rest consumers)
{
    using Batch = std::vector<std::string>;     // ACTG-only chunks pre-split by producer

    constexpr size_t MAX_QUEUE  = 32;           // max batches in flight
    constexpr size_t BATCH_SEQS = 512;          // sequences per batch

    const size_t n_parts     = cfg.num_partitions;
    const size_t n_consumers = n_threads - 1;

    std::deque<Batch>       queue;
    std::mutex              q_mutex;
    std::condition_variable q_cv;
    bool                    producer_done = false;

    std::vector<std::mutex> bucket_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    // Producer: decompress gz one sequence at a time, split on non-ACTG characters,
    // accumulate ACTG-only chunks into batches and push to the queue.
    auto producer_fn = [&]() {
        SeqReader reader;
        if (!reader.load(gz_path)) {
            std::lock_guard<std::mutex> lk(q_mutex);
            producer_done = true;
            q_cv.notify_all();
            return;
        }

        while (true) {
            Batch batch;
            size_t seq_count = 0;
            while (seq_count < BATCH_SEQS && reader.read_next_seq()) {
                split_actg(reader.seq(), reader.seq_len(),
                           [&](const char* c, size_t cl) { batch.emplace_back(c, cl); });
                ++seq_count;
            }
            if (batch.empty()) break;

            {
                std::unique_lock<std::mutex> lk(q_mutex);
                q_cv.wait(lk, [&]{ return queue.size() < MAX_QUEUE; });
                queue.push_back(std::move(batch));
            }
            q_cv.notify_one();
        }

        {
            std::lock_guard<std::mutex> lk(q_mutex);
            producer_done = true;
        }
        q_cv.notify_all();
    };

    // Consumer: pull batches from the queue and extract superkmers.
    auto consumer_fn = [&]() {
        // Cap per-writer buffer so total writer memory ≤ 64 MB/thread regardless of n_parts.
        const size_t flush_thresh = std::max(size_t(4u << 10), size_t(64u << 20) / n_parts);
        MinimizerWindow<k, l>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
        uint64_t local_seqs = 0, local_kmers = 0;

        while (true) {
            Batch batch;
            {
                std::unique_lock<std::mutex> lk(q_mutex);
                q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                if (queue.empty()) break;       // producer_done && queue empty → done
                batch = std::move(queue.front());
                queue.pop_front();
            }
            q_cv.notify_one();                  // wake producer if it was waiting for space

            for (const auto& chunk : batch) {
                extract_superkmers_from_actg<k, l>(
                    chunk.data(), chunk.size(), partition_fn,
                    min_it, writers, local_kmers);
                ++local_seqs;
                for (size_t p = 0; p < n_parts; ++p)
                    if (writers[p].needs_flush())
                        writers[p].flush_to(buckets[p], bucket_mutexes[p]);
            }
        }

        for (size_t p = 0; p < n_parts; ++p)
            writers[p].flush_to(buckets[p], bucket_mutexes[p]);

        total_seqs .fetch_add(local_seqs,  std::memory_order_relaxed);
        total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();

    return { total_seqs.load(), total_kmers.load() };
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
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);

    // Single gz file + multiple requested threads → use producer-consumer.
    // The standard work-stealing path would cap to min(n_threads,1)=1 thread.
    if (n_files == 1 && n_threads_req > 1) {
        const auto& f = cfg.input_files[0];
        if (f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
            return partition_kmers_gz_pc<k, l>(cfg, f, buckets,
                                               partition_fn, n_threads_req);
    }

    const size_t n_threads = std::min(n_threads_req, n_files);
    const size_t n_parts   = cfg.num_partitions;

    std::atomic<size_t>   next_file{0};
    std::vector<std::mutex> bucket_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    auto worker = [&]() {
        // Cap per-writer buffer so total writer memory ≤ 64 MB/thread regardless of n_parts.
        const size_t flush_thresh = std::max(size_t(4u << 10), size_t(64u << 20) / n_parts);
        SeqSource            source;
        MinimizerWindow<k, l>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
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


// ─── Public ───────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t l>
PartitionStats partition_kmers(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets)
{
    const size_t n = cfg.num_partitions;
    return partition_kmers_impl<k, l>(cfg, buckets,
        [n](uint64_t h) -> size_t { return h % n; });
}
