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
#include <stdexcept>
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
    const uint64_t avail = cfg.ram_budget_bytes > 0
        ? cfg.ram_budget_bytes : available_ram_bytes();
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
        if (!cfg.output_kff && !tsv_out) {
            std::cerr << "tuna: error: failed while writing output file: " << cfg.output_file << "\n";
            return 1;
        }

        const double t_phase2 = elapsed_s(t2);
        if (!cfg.hide_progress)
            std::cerr << "      " << total_inserted << " k-mers in  "
                      << total_written << " unique out  " << fmt_s(t_phase2) << "\n";

        std::cerr << "phase1: "        << t_phase1            << "s\n"
                  << "phase2: "        << t_phase2            << "s\n"
                  << "superkmers: "    << stats.superkmers     << "\n"
                  << "n_parts: "       << cfg.num_partitions   << "\n"
                  << "unique_kmers: "  << total_written        << "\n";
        if (!cfg.hide_progress)
            std::cerr << "total: " << fmt_s(elapsed_s(t_start)) << "\n";
        return 0;
    }

    // ── Disk pipeline (fallback when RAM is insufficient) ──────────────────

    std::vector<std::ofstream> buckets(cfg.num_partitions);
    for (size_t p = 0; p < cfg.num_partitions; ++p) {
        const std::string path = partition_path(cfg.work_dir, p);
        buckets[p].open(path, std::ios::binary);
        if (!buckets[p]) {
            std::cerr << "tuna: error: cannot open partition file for writing: " << path << "\n";
            return 1;
        }
    }

    stats = partition_kmers<k, m>(cfg, buckets, disk_write_budget);
    for (size_t p = 0; p < cfg.num_partitions; ++p) {
        buckets[p].close();
        if (!buckets[p]) {
            std::cerr << "tuna: error: failed while writing partition file: "
                      << partition_path(cfg.work_dir, p) << "\n";
            return 1;
        }
    }

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
    if (!cfg.output_kff && !tsv_out) {
        std::cerr << "tuna: error: failed while writing output file: " << cfg.output_file << "\n";
        return 1;
    }

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
    const uint64_t avail   = cfg.ram_budget_bytes > 0
        ? cfg.ram_budget_bytes : available_ram_bytes();
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
        for (size_t p = 0; p < cfg.num_partitions; ++p) {
            const std::string path = partition_path(cfg.work_dir, p);
            buckets[p].open(path, std::ios::binary);
            if (!buckets[p])
                throw std::runtime_error("tuna: cannot open partition file for writing: " + path);
        }
        stats = partition_kmers<k, m>(cfg, buckets, disk_write_budget);
        for (size_t p = 0; p < cfg.num_partitions; ++p) {
            buckets[p].close();
            if (!buckets[p]) {
                throw std::runtime_error(
                    "tuna: failed while writing partition file: " + partition_path(cfg.work_dir, p));
            }
        }

        count_and_callback<k, m>(cfg, stats.kmers, std::forward<Callback>(cb));
    }
}


// ─── Runtime dispatch (k × m) ────────────────────────────────────────────────
//
// Four dispatch modes, selected at compile time:
//
//   FIXED_K + FIXED_M  (recommended for any non-standard k/m pair)
//     Single template instantiation.  Any k ∈ [1,256], any m ∈ [1,k-1].
//     cmake: -DFIXED_K=<k> -DFIXED_M=<m>
//
//   FIXED_K only
//     Instantiates one k with m ∈ {9,11,...,min(k-2,31)} (odd).
//     cmake: -DFIXED_K=<k>
//
//   TUNA_CONDA_PROFILE  (for conda binary packaging)
//     Curated subset: k ∈ {21,31,51,63,127}, ~5 odd m values each (~27 total).
//     Same zero-overhead dispatch as the default build, ~10× smaller binary.
//     cmake: -DTUNA_CONDA_PROFILE=ON
//
//   No flags (default build)
//     Supports k ∈ {11,13,...,127} (odd), m ∈ {9,...,min(k-2,31)} (odd).
//     For any other k/m, prints a message telling the user to recompile.
//
// dispatch_generic is called by both the CLI and the library API.
// It invokes f.template operator()<K, L>() for the matching (K, L) pair.

template <typename F>
inline int dispatch_generic(uint16_t k, uint16_t m, F&& f)
{
#if defined(FIXED_K) && defined(FIXED_M)
    // ── Mode 1: both k and m fixed at compile time ────────────────────────────
    if (k != FIXED_K) {
        std::cerr << "tuna: error: binary compiled for k=" << FIXED_K
                  << "; received k=" << k << "\n"
                  << "  Recompile with -DFIXED_K=" << k << " -DFIXED_M=" << m << "\n";
        return 1;
    }
    if (m != FIXED_M) {
        std::cerr << "tuna: error: binary compiled for m=" << FIXED_M
                  << "; received m=" << m << "\n"
                  << "  Recompile with -DFIXED_K=" << k << " -DFIXED_M=" << m << "\n";
        return 1;
    }
    f.template operator()<FIXED_K, FIXED_M>();
    return 0;

#elif defined(TUNA_CONDA_PROFILE)
    // ── Conda profile: curated k × m subset ──────────────────────────────────
    // k ∈ {21, 31, 51, 63, 127}; odd m values most used in practice.
    // ~27 instantiations total — roughly 10× smaller binary than the full table.
#define TDG(K, L)  if (k == (K) && m == (L)) { f.template operator()<K, L>(); return 0; }
    TDG(21, 11); TDG(21, 13); TDG(21, 15); TDG(21, 17); TDG(21, 19);
    TDG(31, 15); TDG(31, 17); TDG(31, 19); TDG(31, 21); TDG(31, 23);
    TDG(51, 17); TDG(51, 19); TDG(51, 21); TDG(51, 23); TDG(51, 25);
    TDG(63, 19); TDG(63, 21); TDG(63, 23); TDG(63, 25); TDG(63, 27);
    TDG(127, 19); TDG(127, 21); TDG(127, 23); TDG(127, 25); TDG(127, 27); TDG(127, 29); TDG(127, 31);
#undef TDG
    std::cerr << "tuna: error: unsupported k=" << k << " m=" << m << " in conda build\n"
              << "  Conda build supports k ∈ {21,31,51,63,127} with odd m:\n"
              << "    k=21  m ∈ {11,13,15,17,19}\n"
              << "    k=31  m ∈ {15,17,19,21,23}\n"
              << "    k=51  m ∈ {17,19,21,23,25}\n"
              << "    k=63  m ∈ {19,21,23,25,27}\n"
              << "    k=127 m ∈ {19,21,23,25,27,29,31}\n"
              << "  For other k/m recompile with -DFIXED_K=" << k << " -DFIXED_M=" << m << "\n";
    return 1;

#else
    // ── Mode 2 / 3: full dispatch table ──────────────────────────────────────
    // Each k block is guarded by #if so that only the selected k's templates
    // are instantiated — keeping compile time and binary size minimal with FIXED_K.
#define TDG(K, L)  if (k == (K) && m == (L)) { f.template operator()<K, L>(); return 0; }
#if defined(FIXED_K)
    if (k != FIXED_K) {
        std::cerr << "tuna: error: binary compiled for k=" << FIXED_K
                  << "; received k=" << k << "\n"
                  << "  Recompile with -DFIXED_K=" << k << " -DFIXED_M=" << m << "\n";
        return 1;
    }
#endif
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
// k > 32: Kmer<k> spans 2 × uint64_t words.
// m is capped at min(k-2, 31) for the full dispatch build; use -DFIXED_K for
// larger m values to avoid excessive template instantiations.
#if !defined(FIXED_K) || FIXED_K == 33
    TDG(33,  9); TDG(33, 11); TDG(33, 13); TDG(33, 15);
    TDG(33, 17); TDG(33, 19); TDG(33, 21); TDG(33, 23);
    TDG(33, 25); TDG(33, 27); TDG(33, 29); TDG(33, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 35
    TDG(35,  9); TDG(35, 11); TDG(35, 13); TDG(35, 15);
    TDG(35, 17); TDG(35, 19); TDG(35, 21); TDG(35, 23);
    TDG(35, 25); TDG(35, 27); TDG(35, 29); TDG(35, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 37
    TDG(37,  9); TDG(37, 11); TDG(37, 13); TDG(37, 15);
    TDG(37, 17); TDG(37, 19); TDG(37, 21); TDG(37, 23);
    TDG(37, 25); TDG(37, 27); TDG(37, 29); TDG(37, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 39
    TDG(39,  9); TDG(39, 11); TDG(39, 13); TDG(39, 15);
    TDG(39, 17); TDG(39, 19); TDG(39, 21); TDG(39, 23);
    TDG(39, 25); TDG(39, 27); TDG(39, 29); TDG(39, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 41
    TDG(41,  9); TDG(41, 11); TDG(41, 13); TDG(41, 15);
    TDG(41, 17); TDG(41, 19); TDG(41, 21); TDG(41, 23);
    TDG(41, 25); TDG(41, 27); TDG(41, 29); TDG(41, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 43
    TDG(43,  9); TDG(43, 11); TDG(43, 13); TDG(43, 15);
    TDG(43, 17); TDG(43, 19); TDG(43, 21); TDG(43, 23);
    TDG(43, 25); TDG(43, 27); TDG(43, 29); TDG(43, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 45
    TDG(45,  9); TDG(45, 11); TDG(45, 13); TDG(45, 15);
    TDG(45, 17); TDG(45, 19); TDG(45, 21); TDG(45, 23);
    TDG(45, 25); TDG(45, 27); TDG(45, 29); TDG(45, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 47
    TDG(47,  9); TDG(47, 11); TDG(47, 13); TDG(47, 15);
    TDG(47, 17); TDG(47, 19); TDG(47, 21); TDG(47, 23);
    TDG(47, 25); TDG(47, 27); TDG(47, 29); TDG(47, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 49
    TDG(49,  9); TDG(49, 11); TDG(49, 13); TDG(49, 15);
    TDG(49, 17); TDG(49, 19); TDG(49, 21); TDG(49, 23);
    TDG(49, 25); TDG(49, 27); TDG(49, 29); TDG(49, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 51
    TDG(51,  9); TDG(51, 11); TDG(51, 13); TDG(51, 15);
    TDG(51, 17); TDG(51, 19); TDG(51, 21); TDG(51, 23);
    TDG(51, 25); TDG(51, 27); TDG(51, 29); TDG(51, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 53
    TDG(53,  9); TDG(53, 11); TDG(53, 13); TDG(53, 15);
    TDG(53, 17); TDG(53, 19); TDG(53, 21); TDG(53, 23);
    TDG(53, 25); TDG(53, 27); TDG(53, 29); TDG(53, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 55
    TDG(55,  9); TDG(55, 11); TDG(55, 13); TDG(55, 15);
    TDG(55, 17); TDG(55, 19); TDG(55, 21); TDG(55, 23);
    TDG(55, 25); TDG(55, 27); TDG(55, 29); TDG(55, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 57
    TDG(57,  9); TDG(57, 11); TDG(57, 13); TDG(57, 15);
    TDG(57, 17); TDG(57, 19); TDG(57, 21); TDG(57, 23);
    TDG(57, 25); TDG(57, 27); TDG(57, 29); TDG(57, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 59
    TDG(59,  9); TDG(59, 11); TDG(59, 13); TDG(59, 15);
    TDG(59, 17); TDG(59, 19); TDG(59, 21); TDG(59, 23);
    TDG(59, 25); TDG(59, 27); TDG(59, 29); TDG(59, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 61
    TDG(61,  9); TDG(61, 11); TDG(61, 13); TDG(61, 15);
    TDG(61, 17); TDG(61, 19); TDG(61, 21); TDG(61, 23);
    TDG(61, 25); TDG(61, 27); TDG(61, 29); TDG(61, 31);
#endif
#if !defined(FIXED_K) || FIXED_K == 63
    TDG(63,  9); TDG(63, 11); TDG(63, 13); TDG(63, 15);
    TDG(63, 17); TDG(63, 19); TDG(63, 21); TDG(63, 23);
    TDG(63, 25); TDG(63, 27); TDG(63, 29); TDG(63, 31);
#endif
// k=127: Kmer<127> uses 4 × uint64_t words; fully supported by the generic Kmer<k> template.
// m is capped at 31 to avoid excessive instantiations; use -DFIXED_K=127 -DFIXED_M=<m> for larger m.
#if !defined(FIXED_K) || FIXED_K == 127
    TDG(127,  9); TDG(127, 11); TDG(127, 13); TDG(127, 15);
    TDG(127, 17); TDG(127, 19); TDG(127, 21); TDG(127, 23);
    TDG(127, 25); TDG(127, 27); TDG(127, 29); TDG(127, 31);
#endif
#undef TDG

    // k/m pair not in the built-in table.
    std::cerr << "tuna: error: unsupported combination k=" << k << " m=" << m << "\n"
              << "  Default build supports odd k ∈ [11,127] with odd m ∈ [9,min(k-2,31)].\n"
              << "  For any other k or m recompile with:\n"
              << "    -DFIXED_K=" << k << " -DFIXED_M=" << m << "\n";
    return 1;
#endif  // FIXED_K && FIXED_M / TUNA_CONDA_PROFILE / default
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
