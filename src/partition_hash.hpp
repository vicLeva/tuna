#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// Partitioning: minimizer_hash % num_partitions.
//
// Parsing is delegated to SeqSource / GzInput + helicase SIMD parsers, which deliver
// ACTG-only chunks regardless of file format or compression.
//
// Internals:
//   extract_superkmers_from_actg<k, m, PartitionFn>  — pure ACTG sequence logic,
//                                                       no I/O, no threads.
//   partition_kmers_impl<k, m, PartitionFn>           — parallel harness.

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
// minimizer m-mer is tracked as a side-effect of the sliding-window minimum with
// no extra scan.  This replaces the old O(sk_len) find_minimizer_pos rescan that
// consumed ~22% of phase1 cycles.

// FlushFn signature: void(std::vector<SuperkmerWriter>&, size_t partition_id)
// Called immediately after each append — O(1) per superkmer, O(n_superkmers) total.
// Replaces the old O(n_parts) scan loop that ran after every sequence.
template <uint16_t k, uint16_t m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg(
    const char* const             seq,
    const size_t                  seq_len,
    PartitionFn&&                 partition_fn,
    MinimizerWindow<k, m>&        min_it,
    std::vector<SuperkmerWriter>& writers,
    uint64_t&                     kmer_count,
    FlushFn&&                     flush_fn)
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
            flush_fn(writers, pid);
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
    flush_fn(writers, pid);
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

template <uint16_t k, uint16_t m, typename PartitionFn>
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

    // Producer: decompress gz via GzInput, use helicase SIMD parser to deliver
    // ACTG-only chunks, accumulate into batches and push to the queue.
    auto producer_fn = [&]() {
        auto feed = [&](auto& parser) {
            while (true) {
                Batch batch;
                size_t chunk_count = 0;
                while (chunk_count < BATCH_SEQS && parser.next()) {
                    auto [ptr, len] = parser.get_dna_raw();
                    batch.emplace_back(ptr, len);
                    ++chunk_count;
                }
                if (batch.empty()) break;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{ return queue.size() < MAX_QUEUE; });
                    queue.push_back(std::move(batch));
                }
                q_cv.notify_one();
            }
            { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
            q_cv.notify_all();
        };
        try {
            GzInput inp(gz_path);
            if (inp.first_byte() == '@') {
                helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            } else {
                helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(q_mutex);
            producer_done = true;
            q_cv.notify_all();
        }
    };

    // Consumer: pull batches from the queue and extract superkmers.
    auto consumer_fn = [&]() {
        // Cap per-writer buffer so total writer memory ≤ 64 MB/thread regardless of n_parts.
        const size_t flush_thresh = std::max(size_t(4u << 10), size_t(64u << 20) / n_parts);
        MinimizerWindow<k, m>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
        uint64_t local_seqs = 0, local_kmers = 0;

        auto flush_fn = [&](std::vector<SuperkmerWriter>& ws, size_t p) {
            if (ws[p].needs_flush()) ws[p].flush_to(buckets[p], bucket_mutexes[p]);
        };

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
                extract_superkmers_from_actg<k, m>(
                    chunk.data(), chunk.size(), partition_fn,
                    min_it, writers, local_kmers, flush_fn);
                ++local_seqs;
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

template <uint16_t k, uint16_t m, typename PartitionFn>
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
            return partition_kmers_gz_pc<k, m>(cfg, f, buckets,
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
        MinimizerWindow<k, m>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
        uint64_t local_seqs = 0, local_kmers = 0;

        auto flush_fn = [&](std::vector<SuperkmerWriter>& ws, size_t p) {
            if (ws[p].needs_flush()) ws[p].flush_to(buckets[p], bucket_mutexes[p]);
        };

        while (true) {
            const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
            if (fi >= n_files) break;

            source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                extract_superkmers_from_actg<k, m>(
                    chunk, len, partition_fn, min_it, writers, local_kmers, flush_fn);
                ++local_seqs;
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


// ─── Parallel harness (in-memory sinks) ───────────────────────────────────────
//
// Same logic as partition_kmers_impl but writes packed superkmers into
// per-partition std::string buffers instead of disk files.  Used by the
// streaming pipeline to avoid the Phase 1 disk write + Phase 2 mmap round-trip.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_impl(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts      = cfg.num_partitions;

    // Single gz file + multiple threads → producer-consumer (reuse the
    // same consumers but flush to bufs instead of ofstreams).
    // For simplicity we use the multi-file path capped at n_files threads.
    const size_t n_threads = (n_files == 1 && n_threads_req > 1)
        ? n_threads_req   // gz single-file: all threads participate
        : std::min(n_threads_req, n_files);

    // Single .gz file: producer-consumer variant using in-memory sinks.
    if (n_files == 1 && n_threads_req > 1) {
        const auto& gz_path = cfg.input_files[0];
        const bool is_gz = gz_path.size() > 3 &&
                           gz_path.compare(gz_path.size() - 3, 3, ".gz") == 0;
        if (is_gz) {
            using Batch = std::vector<std::string>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t BATCH_SEQS = 512;
            const size_t n_consumers = n_threads - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

            auto producer_fn = [&]() {
                auto feed = [&](auto& parser) {
                    while (true) {
                        Batch batch;
                        size_t chunk_count = 0;
                        while (chunk_count < BATCH_SEQS && parser.next()) {
                            auto [ptr, len] = parser.get_dna_raw();
                            batch.emplace_back(ptr, len);
                            ++chunk_count;
                        }
                        if (batch.empty()) break;
                        { std::unique_lock<std::mutex> lk(q_mutex);
                          q_cv.wait(lk, [&]{ return queue.size() < MAX_QUEUE; });
                          queue.push_back(std::move(batch)); }
                        q_cv.notify_one();
                    }
                    { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
                    q_cv.notify_all();
                };
                try {
                    GzInput inp(gz_path);
                    if (inp.first_byte() == '@') {
                        helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                        feed(p);
                    } else {
                        helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                        feed(p);
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    producer_done = true;
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                const size_t flush_thresh = std::max(size_t(4u << 10), size_t(64u << 20) / n_parts);
                MinimizerWindow<k, m>        min_it;
                std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
                uint64_t local_seqs = 0, local_kmers = 0;
                auto flush_fn = [&](std::vector<SuperkmerWriter>& ws, size_t p) {
                    if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                };
                while (true) {
                    Batch batch;
                    { std::unique_lock<std::mutex> lk(q_mutex);
                      q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                      if (queue.empty()) break;
                      batch = std::move(queue.front());
                      queue.pop_front(); }
                    q_cv.notify_one();
                    for (const auto& chunk : batch) {
                        extract_superkmers_from_actg<k, m>(
                            chunk.data(), chunk.size(), partition_fn,
                            min_it, writers, local_kmers, flush_fn);
                        ++local_seqs;
                    }
                }
                for (size_t p = 0; p < n_parts; ++p)
                    writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
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
    }

    // Multi-file (or single plain file): file-level work-stealing.
    std::atomic<size_t>     next_file{0};
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    auto worker = [&]() {
        const size_t flush_thresh = std::max(size_t(4u << 10), size_t(64u << 20) / n_parts);
        SeqSource            source;
        MinimizerWindow<k, m>        min_it;
        std::vector<SuperkmerWriter> writers(n_parts, SuperkmerWriter(flush_thresh));
        uint64_t local_seqs = 0, local_kmers = 0;

        auto flush_fn = [&](std::vector<SuperkmerWriter>& ws, size_t p) {
            if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
        };

        while (true) {
            const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
            if (fi >= n_files) break;

            source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                extract_superkmers_from_actg<k, m>(
                    chunk, len, partition_fn, min_it, writers, local_kmers, flush_fn);
                ++local_seqs;
            });
        }

        for (size_t p = 0; p < n_parts; ++p)
            writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

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

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets)
{
    const size_t n = cfg.num_partitions;
    return partition_kmers_impl<k, m>(cfg, buckets,
        [n](uint64_t h) -> size_t { return h % n; });
}

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_mem(
    const Config&             cfg,
    std::vector<std::string>& bufs)
{
    const size_t n = cfg.num_partitions;
    return partition_kmers_mem_impl<k, m>(cfg, bufs,
        [n](uint64_t h) -> size_t { return h % n; });
}
