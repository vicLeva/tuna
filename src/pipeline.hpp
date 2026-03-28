#pragma once

// Orchestration only — wire Phase 1 and Phase 2+3 together.
//
// To change what each phase does, edit the corresponding brick header:
//   partition_hash.hpp  — Phase 1: sequences → per-partition superkmer files
//   count.hpp           — Phase 2+3: count k-mers + write output
//   superkmer_io.hpp    — on-disk partition file format

#include "Config.hpp"
#include "kff_output.hpp"
#include "partition_hash.hpp"
#include "count.hpp"

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#ifdef __linux__
#  include <sys/sysinfo.h>
#endif

// Returns available RAM in bytes (Linux sysinfo; 0 if unavailable).
inline uint64_t available_ram_bytes()
{
#ifdef __linux__
    struct sysinfo si{};
    if (sysinfo(&si) == 0)
        return static_cast<uint64_t>(si.freeram + si.bufferram) * si.mem_unit;
#endif
    return 0;
}


// Estimated expansion factor for gzipped genomic FASTA/FASTQ.
// Used to approximate uncompressed size from compressed file size when
// auto-tuning partition count and deciding between in-memory vs disk pipelines.
constexpr uint64_t GZ_EXPAND = 6;

// ─── Timing helper ────────────────────────────────────────────────────────────

inline double elapsed_s(std::chrono::steady_clock::time_point t0)
{
    using Sec = std::chrono::duration<double>;
    return Sec(std::chrono::steady_clock::now() - t0).count();
}

// Format a duration as "X.XXXs" with 3 decimal places.
inline std::string fmt_s(double s)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3fs", s);
    return buf;
}


// ─── Core run function ────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
int run(const Config& cfg)
{
    const auto t_start = std::chrono::steady_clock::now();

    // ── Decide pipeline mode: in-memory vs disk ───────────────────────────
    //
    // In-memory: Phase 1 writes superkmers to per-partition std::string buffers;
    // Phase 2 reads via MemoryReader — no disk round-trip.
    // Selected when estimated packed superkmer size < 60% of available RAM.
    //
    // Estimation: gz → fsz × GZ_EXPAND × 35%; plain → fsz × 2
    // (superkmers store k-1 overlapping bases per boundary ≈ 2× raw sequence).
    uint64_t est_packed = 0;
    for (const auto& f : cfg.input_files) {
        std::error_code ec;
        uint64_t fsz = std::filesystem::file_size(f, ec);
        if (ec) continue;
        const bool is_gz = f.size() > 3 && f.compare(f.size()-3, 3, ".gz") == 0;
        est_packed += is_gz
            ? static_cast<uint64_t>(fsz * GZ_EXPAND * 35 / 100)
            : fsz * 2;
    }
    const uint64_t avail      = available_ram_bytes();
    const bool use_mem_pipeline = avail > 0 && est_packed < avail * 6 / 10;

    // Disk-mode write-buffer budget: up to 30% of RAM across all threads,
    // capped at 512 MB/thread.  Larger buffers → fewer, bigger writes.
    // Falls back to 64 MB/thread if RAM is not detectable.
    const size_t disk_write_budget = [&]() -> size_t {
        if (avail == 0) return size_t(64) << 20;
        const uint64_t budget = avail * 3 / 10 / cfg.num_threads;
        return static_cast<size_t>(std::min(budget, uint64_t(512) << 20));
    }();

    // ── Phase 1: partition ─────────────────────────────────────────────────

    const size_t p1_threads = static_cast<size_t>(cfg.num_threads);
    if (!cfg.hide_progress)
        std::cerr << "[1/2] partitioning  (" << p1_threads << " thread"
                  << (p1_threads > 1 ? "s" : "")
                  << (use_mem_pipeline ? ", in-memory" : "") << ") ...\n";

    PartitionStats stats;
    const auto t_part = std::chrono::steady_clock::now();

    if (use_mem_pipeline) {
        // Pre-reserve buffers to estimated size (×1.2 slack).
        // reserve() is demand-zero — pages are faulted in on first write.
        std::vector<std::string> part_bufs(cfg.num_partitions);
        // Pre-reserve only when per-partition buffers are large enough to matter.
        if (est_packed > 0) {
            const size_t per_part = static_cast<size_t>(
                est_packed * 6 / 5 / cfg.num_partitions);   // ×1.2 slack
            if (per_part >= 64) {
                for (size_t p = 0; p < cfg.num_partitions; ++p)
                    part_bufs[p].reserve(per_part);
            }
        }
        stats = partition_kmers_mem<k, m>(cfg, part_bufs);

        const double t_phase1 = elapsed_s(t_part);
        if (!cfg.hide_progress)
            std::cerr << "      " << stats.seqs << " seqs  "
                      << stats.kmers << " k-mers  " << fmt_s(t_phase1) << "\n";

        if (cfg.partition_only) {
            std::cerr << "phase1: "     << t_phase1        << "s\n"
                      << "superkmers: " << stats.superkmers << "\n";
            if (!cfg.hide_progress)
                std::cerr << "done  (partition only)\n";
            return 0;
        }

        const size_t p2_threads = std::min(static_cast<size_t>(cfg.num_threads),
                                           static_cast<size_t>(cfg.num_partitions));
        if (!cfg.hide_progress)
            std::cerr << "[2/2] counting  (" << p2_threads << " thread"
                      << (p2_threads > 1 ? "s" : "") << ", in-memory) ...\n";

        const auto t2 = std::chrono::steady_clock::now();

        std::ofstream tsv_out;
        if (!cfg.output_kff) {
            tsv_out.open(cfg.output_file);
            if (!tsv_out) {
                std::cerr << "tuna: error: cannot open output file: " << cfg.output_file << "\n";
                return 1;
            }
        }

        const auto [total_inserted, total_written] = cfg.output_kff
            ? [&]() {
                KffOutput kff_out(cfg.output_file, cfg.k);
                auto r = count_and_write_mem<k, m>(cfg, stats.kmers, part_bufs, nullptr, &kff_out);
                kff_out.close();
                return r;
              }()
            : count_and_write_mem<k, m>(cfg, stats.kmers, part_bufs, &tsv_out, nullptr);

        const double t_phase2 = elapsed_s(t2);
        if (!cfg.hide_progress)
            std::cerr << "      " << total_inserted << " k-mers in  "
                      << total_written << " unique out  " << fmt_s(t_phase2) << "\n";

        std::cerr << "phase1: "     << t_phase1        << "s\n"
                  << "phase2: "     << t_phase2        << "s\n"
                  << "superkmers: " << stats.superkmers << "\n";
        if (!cfg.hide_progress)
            std::cerr << "total: " << fmt_s(elapsed_s(t_start)) << "\n";
        return 0;
    }

    // ── Disk pipeline (fallback when RAM is insufficient) ──────────────────

    std::vector<std::ofstream> buckets(cfg.num_partitions);
    for (size_t p = 0; p < cfg.num_partitions; ++p)
        buckets[p].open(partition_path(cfg.work_dir, p),
                        std::ios::binary);

    stats = partition_kmers<k, m>(cfg, buckets, disk_write_budget);
    for (auto& f : buckets) f.close();

    const double t_phase1 = elapsed_s(t_part);
    if (!cfg.hide_progress)
        std::cerr << "      " << stats.seqs << " seqs  "
                  << stats.kmers << " k-mers  " << fmt_s(t_phase1) << "\n";

    if (cfg.partition_only) {
        std::cerr << "phase1: " << t_phase1 << "s\n";
        if (!cfg.hide_progress)
            std::cerr << "done  (partition only)\n";
        return 0;
    }

    // ── Phase 2+3: count + write ───────────────────────────────────────────

    const size_t p2_threads = std::min(static_cast<size_t>(cfg.num_threads),
                                       static_cast<size_t>(cfg.num_partitions));
    if (!cfg.hide_progress)
        std::cerr << "[2/2] counting  (" << p2_threads << " thread"
                  << (p2_threads > 1 ? "s" : "") << ") ...\n";

    const auto t2 = std::chrono::steady_clock::now();

    std::ofstream tsv_out;
    if (!cfg.output_kff) {
        tsv_out.open(cfg.output_file);
        if (!tsv_out) {
            std::cerr << "tuna: error: cannot open output file: " << cfg.output_file << "\n";
            return 1;
        }
    }

    const auto [total_inserted, total_written] = cfg.output_kff
        ? [&]() {
            KffOutput kff_out(cfg.output_file, cfg.k);
            auto r = count_and_write<k, m>(cfg, stats.kmers, nullptr, &kff_out);
            kff_out.close();
            return r;
          }()
        : count_and_write<k, m>(cfg, stats.kmers, &tsv_out, nullptr);

    const double t_phase2 = elapsed_s(t2);
    if (!cfg.hide_progress)
        std::cerr << "      " << total_inserted << " k-mers in  "
                  << total_written << " unique out  " << fmt_s(t_phase2) << "\n";

    // Always emit structured phase timings for tooling (bench.sh).
    std::cerr << "phase1: "     << t_phase1        << "s\n"
              << "phase2: "     << t_phase2        << "s\n"
              << "superkmers: " << stats.superkmers << "\n";

    if (!cfg.hide_progress)
        std::cerr << "total: " << fmt_s(elapsed_s(t_start)) << "\n";

    return 0;
}


// ─── Callback pipeline ────────────────────────────────────────────────────────
//
// run_callback<k, m> mirrors run<k, m> but delivers counted k-mers via a
// user-supplied callback instead of writing to a file.  The same pipeline-mode
// selection (in-memory vs disk) applies.

template <uint16_t k, uint16_t m, typename Callback>
void run_callback(const Config& cfg, Callback&& cb)
{
    // ── Decide pipeline mode (same logic as run<k, m>) ────────────────────
    uint64_t est_packed = 0;
    for (const auto& f : cfg.input_files) {
        std::error_code ec;
        uint64_t fsz = std::filesystem::file_size(f, ec);
        if (ec) continue;
        const bool is_gz = f.size() > 3 && f.compare(f.size()-3, 3, ".gz") == 0;
        est_packed += is_gz
            ? static_cast<uint64_t>(fsz * GZ_EXPAND * 35 / 100)
            : fsz * 2;
    }
    const uint64_t avail   = available_ram_bytes();
    const bool     use_mem = avail > 0 && est_packed < avail * 6 / 10;

    PartitionStats stats;

    if (use_mem) {
        std::vector<std::string> part_bufs(cfg.num_partitions);
        if (est_packed > 0) {
            const size_t per_part = static_cast<size_t>(
                est_packed * 6 / 5 / cfg.num_partitions);
            if (per_part >= 64)
                for (size_t p = 0; p < cfg.num_partitions; ++p)
                    part_bufs[p].reserve(per_part);
        }
        stats = partition_kmers_mem<k, m>(cfg, part_bufs);
        count_and_callback_mem<k, m>(cfg, stats.kmers, part_bufs,
                                     std::forward<Callback>(cb));
    } else {
        const size_t disk_write_budget = [&]() -> size_t {
            if (avail == 0) return size_t(64) << 20;
            const uint64_t budget = avail * 3 / 10 / cfg.num_threads;
            return static_cast<size_t>(std::min(budget, uint64_t(512) << 20));
        }();

        std::vector<std::ofstream> buckets(cfg.num_partitions);
        for (size_t p = 0; p < cfg.num_partitions; ++p)
            buckets[p].open(partition_path(cfg.work_dir, p), std::ios::binary);
        stats = partition_kmers<k, m>(cfg, buckets, disk_write_budget);
        for (auto& f : buckets) f.close();

        count_and_callback<k, m>(cfg, stats.kmers, std::forward<Callback>(cb));
    }
}


// ─── Runtime dispatch (k × m) ────────────────────────────────────────────────
// kache-hash Kmer<k> encodes bases as 2 bits in a single uint64_t → k ≤ 32.
// Only odd k and odd m (even values admit reverse-complement palindromes).
// Supported k: 11..31 (odd)
// Supported m:  9..k-2 (odd, strictly less than k)
//
// dispatch_generic is a single dispatch table used by both the CLI (dispatch)
// and the library API (dispatch_callback), avoiding duplication.
// It invokes f.template operator()<K, L>() for the matching (K, L) pair.
// Requires C++20 template lambdas at the call site.

template <typename F>
inline int dispatch_generic(uint16_t k, uint16_t m, F&& f)
{
#if defined(FIXED_K)
    if (k != FIXED_K) {
        std::cerr << "tuna: error: binary compiled for k=" << FIXED_K
                  << " only; pass -k " << FIXED_K
                  << " or recompile without -DFIXED_K to support other k values\n";
        return 1;
    }
#endif
    // Each k block is guarded by #if so that only the selected k's templates
    // are instantiated — keeping compile time and binary size minimal with FIXED_K.
#define TDG(K, L)  if (k == (K) && m == (L)) { f.template operator()<K, L>(); return 0; }
#if !defined(FIXED_K) || FIXED_K == 11
    TDG(11,  9);
#endif
#if !defined(FIXED_K) || FIXED_K == 13
    TDG(13,  9); TDG(13, 11);
#endif
#if !defined(FIXED_K) || FIXED_K == 15
    TDG(15,  9); TDG(15, 11); TDG(15, 13);
#endif
#if !defined(FIXED_K) || FIXED_K == 17
    TDG(17,  9); TDG(17, 11); TDG(17, 13); TDG(17, 15);
#endif
#if !defined(FIXED_K) || FIXED_K == 19
    TDG(19,  9); TDG(19, 11); TDG(19, 13); TDG(19, 15);
    TDG(19, 17);
#endif
#if !defined(FIXED_K) || FIXED_K == 21
    TDG(21,  9); TDG(21, 11); TDG(21, 13); TDG(21, 15);
    TDG(21, 17); TDG(21, 19);
#endif
#if !defined(FIXED_K) || FIXED_K == 23
    TDG(23,  9); TDG(23, 11); TDG(23, 13); TDG(23, 15);
    TDG(23, 17); TDG(23, 19); TDG(23, 21);
#endif
#if !defined(FIXED_K) || FIXED_K == 25
    TDG(25,  9); TDG(25, 11); TDG(25, 13); TDG(25, 15);
    TDG(25, 17); TDG(25, 19); TDG(25, 21); TDG(25, 23);
#endif
#if !defined(FIXED_K) || FIXED_K == 27
    TDG(27,  9); TDG(27, 11); TDG(27, 13); TDG(27, 15);
    TDG(27, 17); TDG(27, 19); TDG(27, 21); TDG(27, 23);
    TDG(27, 25);
#endif
#if !defined(FIXED_K) || FIXED_K == 29
    TDG(29,  9); TDG(29, 11); TDG(29, 13); TDG(29, 15);
    TDG(29, 17); TDG(29, 19); TDG(29, 21); TDG(29, 23);
    TDG(29, 25); TDG(29, 27);
#endif
#if !defined(FIXED_K) || FIXED_K == 31
    TDG(31,  9); TDG(31, 11); TDG(31, 13); TDG(31, 15);
    TDG(31, 17); TDG(31, 19); TDG(31, 21); TDG(31, 23);
    TDG(31, 25); TDG(31, 27); TDG(31, 29);
#endif
#undef TDG

    std::cerr << "tuna: error: unsupported combination k=" << k << " m=" << m << "\n"
              << "  k must be odd in [11,31]  (kache-hash Kmer fits in 64 bits, k ≤ 32)\n"
              << "  m must be odd in [9,k-2]  (strictly less than k)\n";
    return 1;
}


// CLI dispatcher: run the full pipeline and return its exit code.
inline int dispatch(uint16_t k, uint16_t m, const Config& cfg)
{
    int rc = 1;
    dispatch_generic(k, m, [&]<uint16_t K, uint16_t L>() { rc = run<K, L>(cfg); });
    return rc;
}

// API dispatcher: run the callback pipeline for the given (k, m) pair.
template <typename Callback>
void dispatch_callback(uint16_t k, uint16_t m, const Config& cfg, Callback&& cb)
{
    dispatch_generic(k, m, [&]<uint16_t K, uint16_t L>() {
        run_callback<K, L>(cfg, std::forward<Callback>(cb));
    });
}
