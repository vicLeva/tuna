#pragma once

// Phase 1 — protein sequences → per-partition superkmer files / in-memory buffers.
//
// Partitioning: minimizer_hash % num_partitions.
// Splits on minimizer hash change; each superkmer carries one minimizer.
// min_pos stored as offset from superkmer start → Phase 2 can skip minimizer recompute.

#include "config_prot.hpp"
#include "superkmer_prot.hpp"
#include "minimizer_prot.hpp"
#include "seq_prot.hpp"

#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace prot {

struct PartitionStats {
    uint64_t total_seqs      = 0;
    uint64_t total_kmers     = 0;
    uint64_t total_superkmers = 0;
};

inline size_t prot_writer_flush_threshold(size_t n_parts, size_t budget_per_thread) {
    constexpr size_t MIN_FLUSH = 4u << 10;
    return std::max(MIN_FLUSH, budget_per_thread / n_parts);
}


// ── Core superkmer extraction (single sequence, no I/O) ───────────────────────
//
// Walks one AA-only sequence, detects superkmer boundaries on minimizer hash
// changes, appends each superkmer to the corresponding writer.

template <uint16_t k, uint16_t m, typename FlushFn>
void extract_superkmers_from_prot(
    const char* const                        seq,
    const size_t                             seq_len,
    size_t                                   n_parts,
    MinimizerWindowProt<k, m>&               win,
    std::vector<ProtSuperKmerWriter<k, m>>&  writers,
    uint64_t&                                kmer_count,
    uint64_t&                                sk_count,
    FlushFn&&                                flush_fn,
    std::vector<uint8_t>&                    enc_buf)
{
    using hdr_t = psk_hdr_t<k, m>;
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());

    if (seq_len < k) return;

    // Pre-encode entire sequence to 5-bit values once.
    enc_buf.resize(seq_len);
    for (size_t i = 0; i < seq_len; ++i)
        enc_buf[i] = encode_aa(seq[i]);

    win.reset(enc_buf.data());
    uint64_t prev_hash    = win.hash();
    uint64_t prev_min_pos = win.min_lmer_pos();
    size_t   pid          = prev_hash & (n_parts - 1);
    size_t   sk_start     = 0;

    for (size_t pos = k; pos < seq_len; ++pos) {
        win.advance(enc_buf[pos]);
        const uint64_t new_hash = win.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            const auto min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
            writers[pid].append(enc_buf.data() + sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            prev_min_pos = win.min_lmer_pos();
            pid          = new_hash & (n_parts - 1);
            sk_start     = pos - (k - 1);
        }
    }

    // Flush the final superkmer.
    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    const auto min_pos = static_cast<hdr_t>(prev_min_pos - sk_start);
    writers[pid].append(enc_buf.data() + sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}


// ── In-memory Phase 1 ─────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_mem(
    const ProtConfig&         cfg,
    std::vector<std::string>& part_bufs)   // one string per partition (pre-sized)
{
    const size_t n_parts   = part_bufs.size();
    const size_t n_threads = cfg.num_threads;
    const size_t n_files   = cfg.input_files.size();

    std::vector<std::mutex> buf_mtxs(n_parts);
    std::atomic<size_t>     file_idx{0};
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

    const size_t flush_thresh =
        prot_writer_flush_threshold(n_parts, size_t(64) << 20);

    std::exception_ptr worker_ex;
    std::mutex         ex_mtx;

    auto worker = [&]() {
        try {
            std::vector<ProtSuperKmerWriter<k, m>> writers(n_parts,
                ProtSuperKmerWriter<k, m>(flush_thresh));
            MinimizerWindowProt<k, m> win;
            std::vector<uint8_t> enc_buf;
            ProtSeqSource src;
            uint64_t seqs = 0, kmers = 0, sks = 0;

            auto flush_fn = [&](auto& ws, size_t pid) {
                if (ws[pid].needs_flush())
                    ws[pid].flush_to_mem(part_bufs[pid], buf_mtxs[pid]);
            };

            for (;;) {
                const size_t fi = file_idx.fetch_add(1);
                if (fi >= n_files) break;
                src.process(cfg.input_files[fi],
                    [&](const char* seq, size_t len) {
                        ++seqs;
                        extract_superkmers_from_prot<k, m>(
                            seq, len, n_parts, win, writers,
                            kmers, sks, flush_fn, enc_buf);
                    });
            }

            // Final flush all writers to partition buffers.
            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(part_bufs[p], buf_mtxs[p]);

            total_seqs       += seqs;
            total_kmers      += kmers;
            total_superkmers += sks;
        } catch (...) {
            std::lock_guard<std::mutex> g(ex_mtx);
            if (!worker_ex) worker_ex = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t i = 0; i < n_threads; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    if (worker_ex) std::rethrow_exception(worker_ex);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ── Disk Phase 1 ──────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_disk(
    const ProtConfig&         cfg,
    const std::string&        /* work_dir */,
    std::vector<std::ofstream>& buckets,
    size_t                    write_budget_per_thread)
{
    const size_t n_parts   = buckets.size();
    const size_t n_threads = cfg.num_threads;
    const size_t n_files   = cfg.input_files.size();

    std::vector<std::mutex> bkt_mtxs(n_parts);
    std::atomic<size_t>     file_idx{0};
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

    const size_t flush_thresh =
        prot_writer_flush_threshold(n_parts, write_budget_per_thread);

    std::exception_ptr worker_ex;
    std::mutex         ex_mtx;

    auto worker = [&]() {
        try {
            std::vector<ProtSuperKmerWriter<k, m>> writers(n_parts,
                ProtSuperKmerWriter<k, m>(flush_thresh));
            MinimizerWindowProt<k, m> win;
            std::vector<uint8_t> enc_buf;
            ProtSeqSource src;
            uint64_t seqs = 0, kmers = 0, sks = 0;

            auto flush_fn = [&](auto& ws, size_t pid) {
                if (ws[pid].needs_flush())
                    ws[pid].flush_to(buckets[pid], bkt_mtxs[pid]);
            };

            for (;;) {
                const size_t fi = file_idx.fetch_add(1);
                if (fi >= n_files) break;
                src.process(cfg.input_files[fi],
                    [&](const char* seq, size_t len) {
                        ++seqs;
                        extract_superkmers_from_prot<k, m>(
                            seq, len, n_parts, win, writers,
                            kmers, sks, flush_fn, enc_buf);
                    });
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to(buckets[p], bkt_mtxs[p]);

            total_seqs       += seqs;
            total_kmers      += kmers;
            total_superkmers += sks;
        } catch (...) {
            std::lock_guard<std::mutex> g(ex_mtx);
            if (!worker_ex) worker_ex = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t i = 0; i < n_threads; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    if (worker_ex) std::rethrow_exception(worker_ex);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

} // namespace prot
