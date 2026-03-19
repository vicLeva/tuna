#pragma once

// RAM mode — single-pass parse + count without writing partition files to disk.
//
// Fan-out approach: the producer broadcasts each ACTG chunk (via shared_ptr,
// zero data copy) to ALL n_threads worker queues.  Worker i scans the full
// chunk for minimizer boundaries but only inserts superkmers where
//   pid % n_threads == i
// into its exclusively-owned tables.  No per-table mutex is needed: each
// table is touched by exactly one thread throughout the insert phase.
//
// Working set per thread: (n_partitions / n_threads) tables — e.g. with t=8,
// n=32: 4 × 64 MB = 256 MB per thread instead of 2 GB shared in the old
// single-queue design.  Eliminates both mutex contention and cache thrashing.
//
// Key optimisation vs. disk mode: MinimizerWindow::hash() already is the
// canonical ntHash of the current superkmer's minimizer, forwarded straight
// to Kmer_Window::init_with_hash_ascii() — saves O(sk_len) work per superkmer.
//
// After all input is consumed, the write phase drains every table to the output
// in parallel (work-stealing across threads).

#include "Config.hpp"
#include "count.hpp"
#include "minimizer_window.hpp"
#include "seq_source.hpp"
#include "nt_hash.hpp"

#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"

#include <vector>
#include <deque>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <string>


// ─── Insert-superkmer brick ───────────────────────────────────────────────────
//
// Insert all k-mers from ASCII superkmer seq[0..len-1] into `table`.
// min_nt_h is the canonical ntHash of the m-minimizer — already available from
// MinimizerWindow::hash().  Caller owns the table exclusively — no locking.

template <uint16_t k, uint16_t m>
void insert_superkmer_ascii(
    const char* seq, uint8_t len, uint64_t min_nt_h,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>& table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token)
{
    if (len < k) return;

    auto inc = [](uint32_t v) { return v + 1; };

    kache_hash::Kmer_Window<k, m> win;
    win.init_with_hash_ascii(seq, min_nt_h);
    table.upsert(win, inc, uint32_t(1), token);

    for (uint8_t i = k; i < len; ++i) {
        win.advance(seq[i]);
        table.upsert(win, inc, uint32_t(1), token);
    }
}


// ─── Sequence-level logic brick (filtered) ────────────────────────────────────
//
// Walk one ACTG-only chunk, detect superkmer boundaries via MinimizerWindow.
// Only inserts superkmers where (pid % n_threads == thread_id).
// Table (pid) is exclusively owned by thread (pid % n_threads) — no mutex.

template <uint16_t k, uint16_t m>
void insert_actg_chunk_filtered(
    const char* seq, size_t seq_len,
    MinimizerWindow<k, m>& min_it,
    std::vector<kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>*>& tables,
    std::vector<typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token>& tokens,
    size_t n_parts, size_t thread_id, size_t n_threads,
    uint64_t& kmer_count)
{
    if (seq_len < k) return;

    min_it.reset(seq);
    uint64_t prev_hash = min_it.hash();
    size_t   pid       = prev_hash % n_parts;
    size_t   sk_start  = 0;
    size_t   sk_end    = k;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance(seq[pos]);
        const uint64_t new_hash = min_it.hash();
        if (new_hash != prev_hash || pos - sk_start >= 255u) {
            if (pid % n_threads == thread_id) {
                const auto sk_len = static_cast<uint8_t>(sk_end - sk_start);
                insert_superkmer_ascii<k, m>(seq + sk_start, sk_len, prev_hash,
                                             *tables[pid], tokens[pid]);
                kmer_count += sk_len - k + 1;
            }
            prev_hash = new_hash;
            pid       = new_hash % n_parts;
            sk_start  = pos - (k - 1);
        }
        sk_end = pos + 1;
    }

    // flush last superkmer
    if (pid % n_threads == thread_id) {
        const auto sk_len = static_cast<uint8_t>(sk_end - sk_start);
        insert_superkmer_ascii<k, m>(seq + sk_start, sk_len, prev_hash,
                                     *tables[pid], tokens[pid]);
        kmer_count += sk_len - k + 1;
    }
}


// ─── Parallel harness (fan-out broadcast) ────────────────────────────────────
//
// One bounded queue per worker.  Producer makes one shared_ptr<Task> per chunk
// and pushes it to every worker queue (shared_ptr copies — no data copy).
// Worker i exclusively owns tables p where (p % n_threads == i): no mutex,
// no atomic on the table itself.
//
// After all workers finish, a second parallel phase writes output from each
// table (work-stealing: threads pick the next unwritten partition atomically).

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> partition_and_count_ram(
    const Config&  cfg,
    std::ofstream& out)
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>;

    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts   = static_cast<size_t>(cfg.num_partitions);

    // ── Allocate tables + tokens ───────────────────────────────────────
    // Table p exclusively owned by thread (p % n_threads) — no mutex needed.
    const size_t init_slots = 1u << 22; // 4M slots per table

    std::vector<std::unique_ptr<table_t>> tables(n_parts);
    for (auto& t : tables)
        t = std::make_unique<table_t>(init_slots, 1);

    std::vector<table_t*> table_ptrs(n_parts);
    for (size_t p = 0; p < n_parts; ++p)
        table_ptrs[p] = tables[p].get();

    // One token per table — accessed exclusively by its owning thread.
    std::vector<typename table_t::Token> tokens(n_parts);
    for (size_t p = 0; p < n_parts; ++p)
        tokens[p] = tables[p]->register_user();

    // ── Per-worker queues ─────────────────────────────────────────────
    static constexpr size_t QUEUE_CAP = 64;
    struct Task { std::vector<char> data; };

    struct WorkerQueue {
        std::mutex              mtx;
        std::condition_variable not_empty, not_full;
        std::deque<std::shared_ptr<Task>> q;
        bool done = false;
    };
    std::vector<WorkerQueue> wqs(n_threads);

    std::atomic<uint64_t> total_kmers{0};

    // ── Workers ───────────────────────────────────────────────────────
    auto worker = [&](size_t tid) {
        MinimizerWindow<k, m> min_it;
        uint64_t local_kmers = 0;
        auto& wq = wqs[tid];

        while (true) {
            std::shared_ptr<Task> task;
            {
                std::unique_lock<std::mutex> lk(wq.mtx);
                wq.not_empty.wait(lk, [&]{ return !wq.q.empty() || wq.done; });
                if (wq.q.empty()) break;
                task = std::move(wq.q.front());
                wq.q.pop_front();
            }
            wq.not_full.notify_one();

            insert_actg_chunk_filtered<k, m>(
                task->data.data(), task->data.size(),
                min_it, table_ptrs, tokens,
                n_parts, tid, n_threads, local_kmers);
        }

        total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);

    // ── Producer: broadcast each chunk to all workers ──────────────────
    SeqSource source;
    for (const auto& path : cfg.input_files) {
        source.process(path, [&](const char* chunk, size_t len) {
            auto task = std::make_shared<Task>();
            task->data.assign(chunk, chunk + len);
            for (size_t t = 0; t < n_threads; ++t) {
                auto& wq = wqs[t];
                std::unique_lock<std::mutex> lk(wq.mtx);
                wq.not_full.wait(lk, [&]{ return wq.q.size() < QUEUE_CAP; });
                wq.q.push_back(task);   // shared_ptr copy — no data copy
                wq.not_empty.notify_one();
            }
        });
    }

    // Signal workers that input is exhausted.
    for (size_t t = 0; t < n_threads; ++t) {
        std::lock_guard<std::mutex> lk(wqs[t].mtx);
        wqs[t].done = true;
    }
    for (size_t t = 0; t < n_threads; ++t)
        wqs[t].not_empty.notify_all();
    for (auto& th : threads) th.join();

    // ── Write phase: parallel over partitions ──────────────────────────
    const size_t w_threads = std::min(n_threads, n_parts);
    std::atomic<size_t>   next_part{0};
    std::atomic<uint64_t> total_written{0};
    std::mutex out_mutex;

    auto writer = [&]() {
        std::string chunk;
        while (true) {
            const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
            if (p >= n_parts) break;
            const uint64_t wrt = write_counts(*tables[p], cfg, chunk, out, out_mutex);
            total_written.fetch_add(wrt, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> wthreads;
    wthreads.reserve(w_threads);
    for (size_t t = 0; t < w_threads; ++t)
        wthreads.emplace_back(writer);
    for (auto& th : wthreads) th.join();

    return { total_kmers.load(), total_written.load() };
}
