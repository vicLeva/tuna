#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// Partitioning: minimizer_hash % num_partitions.
//
// Parsing is delegated to SeqSource / GzInput + helicase SIMD parsers, which deliver
// ACTG-only chunks regardless of file format or compression.
//
// Internals:
//   extract_superkmers_from_actg<k, m, partition_m>  — pure ACTG sequence logic,
//                                                       no I/O, no threads.
//   partition_kmers_impl<k, m, partition_m>           — parallel harness.

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
#include <exception>
#include <algorithm>
#include <chrono>


// Per-writer flush threshold: max(4 KB, budget_per_thread / n_parts).
// Disk mode: budget is RAM-proportional — larger buffers yield fewer, bigger writes.
// In-memory mode: fixed 64 MB budget (writes go to RAM, no I/O benefit).
inline size_t writer_flush_threshold(size_t n_parts, size_t budget_per_thread)
{
    constexpr size_t MIN_FLUSH_BYTES = 4u << 10;  // 4 KB minimum
    return std::max(MIN_FLUSH_BYTES, budget_per_thread / n_parts);
}

// Reusable contiguous batch for the compressed producer/consumer path.
// Sequence data is copied once into an arena and boundaries are stored as
// offsets, avoiding one allocation per parsed ACTG chunk.
class ActgBatch
{
    std::vector<char>   bases_;
    std::vector<size_t> ends_;

public:
    void prepare(size_t chunk_capacity)
    {
        constexpr size_t EXPECTED_BASES_PER_CHUNK = 256;
        bases_.clear();
        ends_.clear();
        if (ends_.capacity() < chunk_capacity)
            ends_.reserve(chunk_capacity);
        const size_t base_capacity = chunk_capacity * EXPECTED_BASES_PER_CHUNK;
        if (bases_.capacity() < base_capacity)
            bases_.reserve(base_capacity);
    }

    void append(const char* data, size_t len)
    {
        bases_.insert(bases_.end(), data, data + len);
        ends_.push_back(bases_.size());
    }

    bool empty() const noexcept { return ends_.empty(); }
    size_t size() const noexcept { return ends_.size(); }
    size_t bases_size() const noexcept { return bases_.size(); }

    std::pair<const char*, size_t> chunk(size_t i) const noexcept
    {
        const size_t begin = i == 0 ? 0 : ends_[i - 1];
        return {bases_.data() + begin, ends_[i] - begin};
    }
};


// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only sequence, detect superkmer boundaries on minimizer hash
// changes, and append each superkmer to the corresponding writer.
// Splits on hash change (not partition change) so every superkmer has one
// partition minimizer.
//
// FlushFn: void(std::vector<SuperkmerWriter<k,m>>&, size_t partition_id)
// Called after each append; O(1) per superkmer.
template <uint16_t k, uint16_t m, uint16_t partition_m,
          typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg(
    const char* const                    seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, partition_m>&     min_it,      // partition minimizer iterator
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn,
    std::vector<uint8_t>&                packed_buf)  // packed sequence plus one sentinel byte
{
    using hdr_t = sk_hdr_t<k, m>;  // superkmer header type (local alias)
    // Max value of hdr_t: flush guard prevents sk_len from ever overflowing hdr_t.
    // In practice 2k-m ≤ hdr_t::max() by construction, so this is a safety net only.
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());

    if (seq_len < k) return;
    kmer_count += seq_len - k + 1;

    // Pack the sequence once so overlapping superkmers can copy byte slices
    // instead of repacking the same bases repeatedly.
    const size_t packed_bytes = (seq_len + 3u) / 4u;
    packed_buf.resize(packed_bytes + 1);
    size_t enc = 0;
    for (; enc + 4 <= seq_len; enc += 4) {
        const uint8_t b0 = ((uint8_t(seq[enc    ]) >> 2) ^ (uint8_t(seq[enc    ]) >> 1)) & 3u;
        const uint8_t b1 = ((uint8_t(seq[enc + 1]) >> 2) ^ (uint8_t(seq[enc + 1]) >> 1)) & 3u;
        const uint8_t b2 = ((uint8_t(seq[enc + 2]) >> 2) ^ (uint8_t(seq[enc + 2]) >> 1)) & 3u;
        const uint8_t b3 = ((uint8_t(seq[enc + 3]) >> 2) ^ (uint8_t(seq[enc + 3]) >> 1)) & 3u;
        packed_buf[enc >> 2] = static_cast<uint8_t>(
            (b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
    }
    if (enc < seq_len) {
        uint8_t tail = 0;
        for (size_t i = enc; i < seq_len; ++i) {
            const uint8_t b = ((uint8_t(seq[i]) >> 2) ^ (uint8_t(seq[i]) >> 1)) & 3u;
            tail |= static_cast<uint8_t>(b << (6u - 2u * (i - enc)));
        }
        packed_buf[enc >> 2] = tail;
    }
    packed_buf[packed_bytes] = 0;

    min_it.reset(seq);   // initialises from ASCII (called once per sequence — cheap)
    uint64_t prev_hash = min_it.hash();
    size_t   pid          = partition_fn(prev_hash);  // partition id
    size_t   sk_start     = 0;                        // superkmer start position in current sequence

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance(seq[pos]);
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            const auto min_pos = sk_no_min<k, m>;
            writers[pid].append_packed(packed_buf.data(), sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            ++sk_count;
            prev_hash    = new_hash;
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    const auto min_pos = sk_no_min<k, m>;
    writers[pid].append_packed(packed_buf.data(), sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    ++sk_count;
}


// ─── Producer-consumer harness for compressed files ──────────────────────────
//
// Parser producers feed a shared queue; remaining threads extract and partition
// superkmers. Compressed FASTA uses synchronous producers because partitioning
// dominates CPU. Other compressed formats overlap each parser with an
// AsyncGzInput decompression thread.
//
// The queue is bounded (MAX_QUEUE batches) for backpressure.  Consumers flush
// SuperkmerWriters to the shared bucket files under per-bucket mutexes.

template <uint16_t k, uint16_t m, uint16_t partition_m,
          typename PartitionFn, typename FlushWriterFn>
PartitionStats partition_kmers_gz_pc_impl(
    const Config&               cfg,
    const std::vector<std::string>& gz_paths,
    PartitionFn                 partition_fn,
    FlushWriterFn               flush_writer,
    size_t                      n_threads,          // ≥ 2 (1 producer + rest consumers)
    size_t                      write_budget_per_thread)
{
    using Batch = ActgBatch;

    constexpr size_t MAX_QUEUE   = 32;       // max batches in flight
    constexpr size_t BATCH_SEQS  = 512;      // sequences per batch
    constexpr size_t BATCH_BASES = 8u << 20; // keep long FASTA batches parallel
    constexpr size_t DYNAMIC_JOIN_BASES = 32u << 20;
    constexpr auto BACKPRESSURE_WAIT = std::chrono::microseconds(50);

    const size_t n_parts = cfg.num_partitions;
    const bool all_gz_fasta = std::all_of(
        gz_paths.begin(), gz_paths.end(),
        [](const std::string& path) {
            return path.ends_with(".fa.gz") ||
                   path.ends_with(".fasta.gz") ||
                   path.ends_with(".fna.gz");
        });
    // Start with up to one input-capable worker per four requested threads.
    // FASTA inflation stays synchronous because partitioning dominates on
    // assemblies. Read files overlap each parser with an inflater, then reclaim
    // that thread-budget slot for counting when its input is exhausted.
    const bool async_input = n_threads >= 3 && !all_gz_fasta;
    const size_t n_producers = std::min(
        gz_paths.size(), std::max(size_t(1), n_threads / 4));
    const size_t n_consumers = n_threads -
        n_producers * (async_input ? size_t(2) : size_t(1));

    std::deque<Batch>       queue;
    std::deque<Batch>       free_batches;
    size_t                  queued_bases = 0;
    std::mutex              q_mutex;
    std::condition_variable q_cv;
    size_t                  producers_done = 0;
    size_t                  async_slots_released = 0;
    size_t                  async_helpers_started = 0;
    std::atomic<size_t>     next_gz{0};
    std::exception_ptr      producer_error = nullptr;
    std::atomic<bool>       stop{false};
    std::exception_ptr      consumer_error = nullptr;
    std::mutex              consumer_error_mutex;

    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

    struct WorkerState {
        MinimizerWindow<k, partition_m> min_it;
        std::vector<SuperkmerWriter<k, m>> writers;
        std::vector<uint8_t> packed_buf;
        uint64_t seqs = 0;
        uint64_t kmers = 0;
        uint64_t superkmers = 0;

        WorkerState(size_t parts, size_t flush_thresh)
            : writers(parts, SuperkmerWriter<k, m>(flush_thresh)) {}
    };

    // Input-capable workers create partitioning state lazily. This avoids tens
    // of thousands of per-partition writer allocations when a small read tail
    // cannot amortize recruiting another consumer.
    auto producer_fn = [&]() {
        bool producer_finished = false;
        std::unique_ptr<WorkerState> worker_state;
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            auto ensure_state = [&]() -> WorkerState& {
                if (!worker_state)
                    worker_state = std::make_unique<WorkerState>(n_parts, flush_thresh);
                return *worker_state;
            };

            auto process_batch = [&](Batch&& batch) {
                auto& state = ensure_state();
                auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                    if (ws[p].needs_flush()) flush_writer(ws[p], p);
                };
                for (size_t i = 0; i < batch.size(); ++i) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    const auto [chunk, chunk_len] = batch.chunk(i);
                    extract_superkmers_from_actg<k, m, partition_m>(
                        chunk, chunk_len, partition_fn,
                        state.min_it, state.writers, state.kmers, state.superkmers,
                        flush_fn, state.packed_buf);
                    ++state.seqs;
                }
                {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    free_batches.push_back(std::move(batch));
                }
                q_cv.notify_one();
            };

            auto consume_one = [&]() -> bool {
                Batch batch;
                {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    if (queue.empty()) return false;
                    queued_bases -= queue.front().bases_size();
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                process_batch(std::move(batch));
                return true;
            };

            auto feed = [&](auto& parser) {
                    while (!stop.load(std::memory_order_relaxed)) {
                        Batch batch;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (!free_batches.empty()) {
                                batch = std::move(free_batches.front());
                                free_batches.pop_front();
                            }
                        }
                        batch.prepare(BATCH_SEQS);
                        size_t chunk_count = 0;
                        while (!stop.load(std::memory_order_relaxed)
                               && chunk_count < BATCH_SEQS
                               && batch.bases_size() < BATCH_BASES
                               && parser.next()) {
                            auto [ptr, len] = parser.get_dna_raw();
                            batch.append(ptr, len);
                            ++chunk_count;
                        }
                        if (stop.load(std::memory_order_relaxed) || batch.empty()) break;

                        while (!stop.load(std::memory_order_relaxed)) {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            if (queue.size() < MAX_QUEUE) {
                                queued_bases += batch.bases_size();
                                queue.push_back(std::move(batch));
                                lk.unlock();
                                q_cv.notify_one();
                                break;
                            }

                            // A sole async parser is input-critical for a
                            // single compressed stream. Keep it dedicated.
                            if (async_input && n_producers == 1) {
                                q_cv.wait(lk, [&] {
                                    return stop.load(std::memory_order_relaxed) ||
                                           queue.size() < MAX_QUEUE;
                                });
                                continue;
                            }

                            // Ignore transient fullness. Sustained backpressure
                            // indicates that this worker is more valuable as a
                            // partition consumer on the current workload.
                            const bool ready = q_cv.wait_for(
                                lk, BACKPRESSURE_WAIT, [&] {
                                    return stop.load(std::memory_order_relaxed) ||
                                           queue.size() < MAX_QUEUE;
                                });
                            if (ready) continue;
                            lk.unlock();
                            if (!consume_one()) std::this_thread::yield();
                        }
                    }
            };

            while (!stop.load(std::memory_order_relaxed)) {
                const size_t file_idx =
                    next_gz.fetch_add(1, std::memory_order_relaxed);
                if (file_idx >= gz_paths.size()) break;
                const auto& gz_path = gz_paths[file_idx];
                if (async_input) {
                    AsyncGzInput inp(gz_path);
                    if (inp.first_byte() == '@') {
                        helicase::FastqParser<HELICASE_ACTG, AsyncGzInput> p(std::move(inp));
                        feed(p);
                    } else {
                        helicase::FastaParser<HELICASE_ACTG, AsyncGzInput> p(std::move(inp));
                        feed(p);
                    }
                } else {
                    GzInput inp(gz_path);
                    if (inp.first_byte() == '@') {
                        helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                        feed(p);
                    } else {
                        helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                        feed(p);
                    }
                }
            }

            {
                std::lock_guard<std::mutex> lk(q_mutex);
                ++producers_done;
                if (async_input) ++async_slots_released;
            }
            producer_finished = true;
            q_cv.notify_all();

            // Drain until every producer has finished and no queued work remains.
            bool should_drain = true;
            if (!worker_state) {
                std::lock_guard<std::mutex> lk(q_mutex);
                should_drain = queued_bases >= DYNAMIC_JOIN_BASES;
            }
            if (should_drain) {
                while (!stop.load(std::memory_order_relaxed)) {
                    Batch batch;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&] {
                            return stop.load(std::memory_order_relaxed) ||
                                   !queue.empty() || producers_done == n_producers;
                        });
                        if (stop.load(std::memory_order_relaxed)) break;
                        if (queue.empty() && producers_done == n_producers) break;
                        if (queue.empty()) continue;
                        queued_bases -= queue.front().bases_size();
                        batch = std::move(queue.front());
                        queue.pop_front();
                    }
                    process_batch(std::move(batch));
                }
            }

            if (worker_state) {
                for (size_t p = 0; p < n_parts; ++p)
                    flush_writer(worker_state->writers[p], p);
                total_seqs.fetch_add(worker_state->seqs, std::memory_order_relaxed);
                total_kmers.fetch_add(worker_state->kmers, std::memory_order_relaxed);
                total_superkmers.fetch_add(
                    worker_state->superkmers, std::memory_order_relaxed);
            }
        } catch (...) {
            if (!producer_finished) {
                std::lock_guard<std::mutex> lk(q_mutex);
                ++producers_done;
                if (async_input) ++async_slots_released;
                producer_finished = true;
            }
            {
                std::lock_guard<std::mutex> lk(consumer_error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    // Ordinary consumers keep their hot state stack-local so the compiler can
    // retain counters and minimizer state in registers.
    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh =
                writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, partition_m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(
                n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> packed_buf;
            uint64_t local_seqs = 0;
            uint64_t local_kmers = 0;
            uint64_t local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) flush_writer(ws[p], p);
            };

            while (!stop.load(std::memory_order_relaxed)) {
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&] {
                        return stop.load(std::memory_order_relaxed) ||
                               !queue.empty() || producers_done == n_producers;
                    });
                    if (stop.load(std::memory_order_relaxed)) break;
                    if (queue.empty() && producers_done == n_producers) break;
                    if (queue.empty()) continue;
                    queued_bases -= queue.front().bases_size();
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();

                for (size_t i = 0; i < batch.size(); ++i) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    const auto [chunk, chunk_len] = batch.chunk(i);
                    extract_superkmers_from_actg<k, m, partition_m>(
                        chunk, chunk_len, partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn,
                        packed_buf);
                    ++local_seqs;
                }
                {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    free_batches.push_back(std::move(batch));
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                flush_writer(writers[p], p);
            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(consumer_error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_producers; ++t)
        threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    // AsyncGzInput consumes one thread-budget slot per producer. Keep its
    // replacement asleep until that producer has finished all input, then use
    // the released slot to help drain the remaining partition backlog.
    if (async_input) {
        for (size_t t = 0; t < n_producers; ++t) {
            threads.emplace_back([&] {
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&] {
                        return stop.load(std::memory_order_relaxed) ||
                               async_helpers_started < async_slots_released;
                    });
                    if (stop.load(std::memory_order_relaxed)) return;
                    ++async_helpers_started;
                    if (queued_bases < DYNAMIC_JOIN_BASES) return;
                }
                consumer_fn();
            });
        }
    }
    for (auto& th : threads) th.join();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, uint16_t partition_m, typename PartitionFn>
PartitionStats partition_kmers_gz_pc(
    const Config&                  cfg,
    const std::vector<std::string>& gz_paths,
    std::vector<std::ofstream>&    buckets,
    PartitionFn                    partition_fn,
    size_t                         n_threads,
    size_t                         write_budget_per_thread)
{
    std::vector<std::mutex> bucket_mutexes(cfg.num_partitions);
    auto flush_writer = [&](SuperkmerWriter<k, m>& writer, size_t p) {
        writer.flush_to(buckets[p], bucket_mutexes[p]);
    };
    return partition_kmers_gz_pc_impl<k, m, partition_m>(
        cfg, gz_paths, partition_fn, flush_writer,
        n_threads, write_budget_per_thread);
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// File-level work stealing: each worker atomically claims the next file and runs
// helicase + superkmer extraction on its own thread.
// Load balancing is coarse (file granularity) but negligible for equal-size files.

template <uint16_t k, uint16_t m, uint16_t partition_m, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);

    const bool all_gz = std::all_of(
        cfg.input_files.begin(), cfg.input_files.end(),
        [](const std::string& f) {
            return f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0;
        });
    if (n_threads_req > 1 && all_gz)
        return partition_kmers_gz_pc<k, m, partition_m>(
            cfg, cfg.input_files, buckets, partition_fn,
            n_threads_req, write_budget_per_thread);

    const size_t n_threads = std::min(n_threads_req, n_files);
    const size_t n_parts   = cfg.num_partitions;

    std::atomic<size_t>   next_file{0};
    std::vector<std::mutex> bucket_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            SeqSource            source;
            MinimizerWindow<k, partition_m>     min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         packed_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to(buckets[p], bucket_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m, partition_m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers,
                        flush_fn, packed_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to(buckets[p], bucket_mutexes[p]);

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Parallel harness (in-memory sinks) ───────────────────────────────────────
//
// Same logic as partition_kmers_impl but writes packed superkmers into
// per-partition std::string buffers instead of disk files.  Used by the
// streaming pipeline to avoid the Phase 1 disk write + Phase 2 mmap round-trip.

template <uint16_t k, uint16_t m, uint16_t partition_m, typename PartitionFn>
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

    const bool all_gz = std::all_of(
        cfg.input_files.begin(), cfg.input_files.end(),
        [](const std::string& f) {
            return f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0;
        });
    if (n_threads_req > 1 && all_gz) {
        std::vector<std::mutex> buf_mutexes(n_parts);
        auto flush_writer = [&](SuperkmerWriter<k, m>& writer, size_t p) {
            writer.flush_to_mem(bufs[p], buf_mutexes[p]);
        };
        return partition_kmers_gz_pc_impl<k, m, partition_m>(
            cfg, cfg.input_files, partition_fn, flush_writer,
            n_threads_req, 64u << 20);
    }

    // Single .gz file: producer-consumer variant using in-memory sinks.
    if (n_files == 1 && n_threads_req > 1) {
        const auto& gz_path = cfg.input_files[0];
        const bool is_gz = gz_path.size() > 3 &&
                           gz_path.compare(gz_path.size() - 3, 3, ".gz") == 0;
        if (is_gz) {
            using Batch = ActgBatch;
            constexpr size_t MAX_QUEUE   = 32;
            constexpr size_t BATCH_SEQS  = 512;
            constexpr size_t BATCH_BASES = 8u << 20;
            const bool async_input = n_threads >= 3;
            const size_t n_consumers = n_threads - (async_input ? 2 : 1);

            std::deque<Batch>       queue;
            std::deque<Batch>       free_batches;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto feed = [&](auto& parser) {
                    while (true) {
                        Batch batch;
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            if (!free_batches.empty()) {
                                batch = std::move(free_batches.front());
                                free_batches.pop_front();
                            }
                        }
                        batch.prepare(BATCH_SEQS);
                        size_t chunk_count = 0;
                        while (!stop.load(std::memory_order_relaxed)
                               && chunk_count < BATCH_SEQS
                               && batch.bases_size() < BATCH_BASES
                               && parser.next()) {
                            auto [ptr, len] = parser.get_dna_raw();
                            batch.append(ptr, len);
                            ++chunk_count;
                        }
                        if (stop.load(std::memory_order_relaxed)) break;
                        if (batch.empty()) break;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&] {
                                return stop.load(std::memory_order_relaxed) ||
                                       queue.size() < MAX_QUEUE;
                            });
                            if (stop.load(std::memory_order_relaxed)) break;
                            queue.push_back(std::move(batch));
                        }
                        q_cv.notify_one();
                    }
                    { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
                    q_cv.notify_all();
                };
                try {
                    if (async_input) {
                        AsyncGzInput inp(gz_path);
                        if (inp.first_byte() == '@') {
                            helicase::FastqParser<HELICASE_ACTG, AsyncGzInput> p(std::move(inp));
                            feed(p);
                        } else {
                            helicase::FastaParser<HELICASE_ACTG, AsyncGzInput> p(std::move(inp));
                            feed(p);
                        }
                    } else {
                        GzInput inp(gz_path);
                        if (inp.first_byte() == '@') {
                            helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                            feed(p);
                        } else {
                            helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                            feed(p);
                        }
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, partition_m> min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         packed_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&] {
                                return stop.load(std::memory_order_relaxed) ||
                                       !queue.empty() || producer_done;
                            });
                            if (stop.load(std::memory_order_relaxed) || queue.empty()) break;
                            batch = std::move(queue.front());
                            queue.pop_front();
                        }
                        q_cv.notify_one();
                        for (size_t i = 0; i < batch.size(); ++i) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            const auto [chunk, chunk_len] = batch.chunk(i);
                            extract_superkmers_from_actg<k, m, partition_m>(
                                chunk, chunk_len, partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn,
                                packed_buf);
                            ++local_seqs;
                        }
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            free_batches.push_back(std::move(batch));
                        }
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
        }
    }

    // Multi-file (or single plain file): file-level work-stealing.
    std::atomic<size_t>     next_file{0};
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            SeqSource            source;
            MinimizerWindow<k, partition_m>     min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         packed_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m, partition_m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers,
                        flush_fn, packed_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);
    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Public ───────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    size_t                      write_budget_per_thread)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    static constexpr uint16_t partition_m = m;
    return partition_kmers_impl<k, m, partition_m>(cfg, buckets,
        [mask](uint64_t h) -> size_t { return h & mask; },
        write_budget_per_thread);
}

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_mem(
    const Config&             cfg,
    std::vector<std::string>& bufs)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    static constexpr uint16_t partition_m = m;
    return partition_kmers_mem_impl<k, m, partition_m>(cfg, bufs,
        [mask](uint64_t h) -> size_t { return h & mask; });
}
