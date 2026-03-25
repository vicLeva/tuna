#pragma once

// Orchestration only — wire Phase 1 and Phase 2+3 together.
//
// To change what each phase does, edit the corresponding brick header:
//   partition_hash.hpp  — Phase 1: sequences → per-partition superkmer files
//   count.hpp           — Phase 2+3: count k-mers + write output
//   superkmer_io.hpp    — on-disk partition file format

#include "Config.hpp"
#include "partition_hash.hpp"
#include "count.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
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


// ─── Timing helper ────────────────────────────────────────────────────────────

inline double elapsed_s(std::chrono::steady_clock::time_point t0)
{
    using Sec = std::chrono::duration<double>;
    return Sec(std::chrono::steady_clock::now() - t0).count();
}


// ─── Core run function ────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
int run(const Config& cfg)
{
    const auto t_start = std::chrono::steady_clock::now();

    // ── Decide pipeline mode: in-memory (no disk round-trip) vs disk ──────
    //
    // In-memory mode: Phase 1 writes packed superkmers to per-partition
    // std::string buffers; Phase 2 reads from them via MemoryReader.
    // Saves the disk write (Phase 1) + mmap setup (Phase 2) overhead.
    //
    // Condition: estimated packed superkmer size < 60% of available RAM.
    // Packed size ≈ 35% of uncompressed sequence size.
    // .gz files are expanded by GZ_EXPAND=6 for the estimate.
    uint64_t est_raw = 0;
    for (const auto& f : cfg.input_files) {
        std::error_code ec;
        uint64_t fsz = std::filesystem::file_size(f, ec);
        if (ec) continue;
        const bool is_gz = f.size() > 3 && f.compare(f.size()-3, 3, ".gz") == 0;
        est_raw += is_gz ? fsz * 6 : fsz;
    }
    const uint64_t est_packed = static_cast<uint64_t>(est_raw * 0.35);
    const uint64_t avail      = available_ram_bytes();
    const bool use_mem_pipeline = avail > 0 && est_packed < avail * 6 / 10;

    // ── Phase 1: partition ─────────────────────────────────────────────────

    const size_t p1_threads = static_cast<size_t>(cfg.num_threads);
    if (!cfg.hide_progress)
        std::cerr << "[1/2] Partitioning into " << cfg.num_partitions
                  << " partitions  (" << p1_threads << " thread"
                  << (p1_threads > 1 ? "s" : "")
                  << (use_mem_pipeline ? ", in-memory" : "") << ")...\n";

    PartitionStats stats;
    const auto t_part = std::chrono::steady_clock::now();

    if (use_mem_pipeline) {
        // In-memory: write to per-partition string buffers.
        // Pre-reserve each buffer to the estimated final size (with 20% slack) and
        // touch all pages upfront — converts 12M+ scattered page faults into a single
        // sequential zero-fill, eliminating the sys-time spike on large partition counts.
        std::vector<std::string> part_bufs(cfg.num_partitions);
        // Pre-reserve only when total partition memory > 2 GB.
        // Below that, buffers fit in L3 and grow cheaply; pre-faulting upfront
        // would merely increase working-set size and add cache pressure.
        // Above that threshold, 32K+ growing strings cause O(millions) scattered
        // page faults (measurable as 100s+ of sys time) that pre-faulting eliminates.
        if (est_packed > 0) {
            const size_t per_part = static_cast<size_t>(
                est_packed * 6 / 5 / cfg.num_partitions);   // ×1.2 slack
            const uint64_t total = static_cast<uint64_t>(per_part) * cfg.num_partitions;
            if (per_part >= 64 && total > (2ULL << 30)) {
                for (auto& s : part_bufs) {
                    s.resize(per_part, '\0');   // allocate + touch all pages
                    s.clear();                   // size → 0, capacity stays
                }
            }
        }
        stats = partition_kmers_mem<k, m>(cfg, part_bufs);

        const double t_phase1 = elapsed_s(t_part);
        if (!cfg.hide_progress)
            std::cerr << "    " << stats.seqs  << " sequences, "
                      << stats.kmers << " k-mers  ("
                      << t_phase1 << "s)\n";

        if (cfg.partition_only) {
            std::cerr << "phase1: "     << t_phase1        << "s\n"
                      << "superkmers: " << stats.superkmers << "\n";
            if (!cfg.hide_progress)
                std::cerr << "[done] partition only\n";
            return 0;
        }

        const size_t p2_threads = std::min(static_cast<size_t>(cfg.num_threads),
                                           static_cast<size_t>(cfg.num_partitions));
        if (!cfg.hide_progress)
            std::cerr << "[2/2] Counting + writing  (" << p2_threads << " thread"
                      << (p2_threads > 1 ? "s" : "") << ", in-memory)...\n";

        const auto t2 = std::chrono::steady_clock::now();
        std::ofstream out(cfg.output_file);
        if (!out) {
            std::cerr << "tuna: error: cannot open output file: " << cfg.output_file << "\n";
            return 1;
        }
        const auto [total_inserted, total_written] =
            count_and_write_mem<k, m>(cfg, stats.kmers, part_bufs, out);

        const double t_phase2 = elapsed_s(t2);
        if (!cfg.hide_progress)
            std::cerr << "    " << total_inserted << " k-mers processed, "
                      << total_written << " unique written  ("
                      << t_phase2 << "s)\n";

        std::cerr << "phase1: "     << t_phase1        << "s\n"
                  << "phase2: "     << t_phase2        << "s\n"
                  << "superkmers: " << stats.superkmers << "\n";
        if (!cfg.hide_progress)
            std::cerr << "[done] total: " << elapsed_s(t_start) << "s\n";
        return 0;
    }

    // ── Disk pipeline (fallback when RAM is insufficient) ──────────────────

    std::vector<std::ofstream> buckets(cfg.num_partitions);
    for (size_t p = 0; p < cfg.num_partitions; ++p)
        buckets[p].open(cfg.work_dir + "hash_" + std::to_string(p) + ".superkmers",
                        std::ios::binary);

    stats = partition_kmers<k, m>(cfg, buckets);
    for (auto& f : buckets) f.close();

    const double t_phase1 = elapsed_s(t_part);
    if (!cfg.hide_progress)
        std::cerr << "    " << stats.seqs  << " sequences, "
                  << stats.kmers << " k-mers  ("
                  << t_phase1 << "s)\n";

    if (cfg.partition_only) {
        std::cerr << "phase1: " << t_phase1 << "s\n";
        if (!cfg.hide_progress)
            std::cerr << "[done] partition only\n";
        return 0;
    }

    // ── Phase 2+3: count + write ───────────────────────────────────────────

    const size_t p2_threads = std::min(static_cast<size_t>(cfg.num_threads),
                                       static_cast<size_t>(cfg.num_partitions));
    if (!cfg.hide_progress)
        std::cerr << "[2/2] Counting + writing  (" << p2_threads << " thread"
                  << (p2_threads > 1 ? "s" : "") << ")...\n";

    const auto t2 = std::chrono::steady_clock::now();

    std::ofstream out(cfg.output_file);
    if (!out) {
        std::cerr << "tuna: error: cannot open output file: " << cfg.output_file << "\n";
        return 1;
    }

    const auto [total_inserted, total_written] = count_and_write<k, m>(cfg, stats.kmers, out);

    const double t_phase2 = elapsed_s(t2);
    if (!cfg.hide_progress)
        std::cerr << "    " << total_inserted << " k-mers processed, "
                  << total_written << " unique written  ("
                  << t_phase2 << "s)\n";

    // Always emit structured phase timings for tooling (bench.sh).
    std::cerr << "phase1: "     << t_phase1        << "s\n"
              << "phase2: "     << t_phase2        << "s\n"
              << "superkmers: " << stats.superkmers << "\n";

    if (!cfg.hide_progress)
        std::cerr << "[done] total: " << elapsed_s(t_start) << "s\n";

    return 0;
}


// ─── Runtime dispatch (k × m) ────────────────────────────────────────────────
// kache-hash Kmer<k> encodes bases as 2 bits in a single uint64_t → k ≤ 32.
// Only odd k and odd m (even values admit reverse-complement palindromes).
// Supported k: 11..31 (odd)
// Supported m:  9..k-2 (odd, strictly less than k)

#define TUNA_DISPATCH(K, L)  if (k == (K) && m == (L)) return run<K, L>(cfg)

inline int dispatch(uint16_t k, uint16_t m, const Config& cfg)
{
    // k = 11  (1 pair)
    TUNA_DISPATCH(11,  9);
    // k = 13  (2 pairs)
    TUNA_DISPATCH(13,  9); TUNA_DISPATCH(13, 11);
    // k = 15  (3 pairs)
    TUNA_DISPATCH(15,  9); TUNA_DISPATCH(15, 11); TUNA_DISPATCH(15, 13);
    // k = 17  (4 pairs)
    TUNA_DISPATCH(17,  9); TUNA_DISPATCH(17, 11); TUNA_DISPATCH(17, 13); TUNA_DISPATCH(17, 15);
    // k = 19  (5 pairs)
    TUNA_DISPATCH(19,  9); TUNA_DISPATCH(19, 11); TUNA_DISPATCH(19, 13); TUNA_DISPATCH(19, 15);
    TUNA_DISPATCH(19, 17);
    // k = 21  (6 pairs)
    TUNA_DISPATCH(21,  9); TUNA_DISPATCH(21, 11); TUNA_DISPATCH(21, 13); TUNA_DISPATCH(21, 15);
    TUNA_DISPATCH(21, 17); TUNA_DISPATCH(21, 19);
    // k = 23  (7 pairs)
    TUNA_DISPATCH(23,  9); TUNA_DISPATCH(23, 11); TUNA_DISPATCH(23, 13); TUNA_DISPATCH(23, 15);
    TUNA_DISPATCH(23, 17); TUNA_DISPATCH(23, 19); TUNA_DISPATCH(23, 21);
    // k = 25  (8 pairs)
    TUNA_DISPATCH(25,  9); TUNA_DISPATCH(25, 11); TUNA_DISPATCH(25, 13); TUNA_DISPATCH(25, 15);
    TUNA_DISPATCH(25, 17); TUNA_DISPATCH(25, 19); TUNA_DISPATCH(25, 21); TUNA_DISPATCH(25, 23);
    // k = 27  (9 pairs)
    TUNA_DISPATCH(27,  9); TUNA_DISPATCH(27, 11); TUNA_DISPATCH(27, 13); TUNA_DISPATCH(27, 15);
    TUNA_DISPATCH(27, 17); TUNA_DISPATCH(27, 19); TUNA_DISPATCH(27, 21); TUNA_DISPATCH(27, 23);
    TUNA_DISPATCH(27, 25);
    // k = 29  (10 pairs)
    TUNA_DISPATCH(29,  9); TUNA_DISPATCH(29, 11); TUNA_DISPATCH(29, 13); TUNA_DISPATCH(29, 15);
    TUNA_DISPATCH(29, 17); TUNA_DISPATCH(29, 19); TUNA_DISPATCH(29, 21); TUNA_DISPATCH(29, 23);
    TUNA_DISPATCH(29, 25); TUNA_DISPATCH(29, 27);
    // k = 31  (11 pairs)
    TUNA_DISPATCH(31,  9); TUNA_DISPATCH(31, 11); TUNA_DISPATCH(31, 13); TUNA_DISPATCH(31, 15);
    TUNA_DISPATCH(31, 17); TUNA_DISPATCH(31, 19); TUNA_DISPATCH(31, 21); TUNA_DISPATCH(31, 23);
    TUNA_DISPATCH(31, 25); TUNA_DISPATCH(31, 27); TUNA_DISPATCH(31, 29);

    std::cerr << "tuna: error: unsupported combination k=" << k << " m=" << m << "\n"
              << "  k must be odd in [11,31]  (kache-hash Kmer fits in 64 bits, k ≤ 32)\n"
              << "  m must be odd in [9,k-2]  (strictly less than k)\n";
    return 1;
}

#undef TUNA_DISPATCH
