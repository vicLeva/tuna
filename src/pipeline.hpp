#pragma once

// Orchestration only — wire Phase 1 and Phase 2+3 together.
//
// To change what each phase does, edit the corresponding brick header:
//   partition_hash.hpp  — hash and kmtricks partition strategies
//   partition_kmc.hpp   — KMC signature partition strategy (default)
//   count.hpp           — how k-mers are counted and output is written
//   superkmer_io.hpp    — the on-disk partition file format

#include "Config.hpp"
#include "partition_hash.hpp"
#include "partition_kmc.hpp"
#include "count.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>


// ─── Timing helper ────────────────────────────────────────────────────────────

inline double elapsed_s(std::chrono::steady_clock::time_point t0)
{
    using Sec = std::chrono::duration<double>;
    return Sec(std::chrono::steady_clock::now() - t0).count();
}


// ─── Core run function ────────────────────────────────────────────────────────

template <uint16_t k, uint16_t l>
int run(const Config& cfg)
{
    const auto t_start = std::chrono::steady_clock::now();

    // ── Phase 1: partition ─────────────────────────────────────────────────

    const size_t p1_threads = std::min(static_cast<size_t>(cfg.num_threads),
                                       cfg.input_files.size());
    if (!cfg.hide_progress)
        std::cerr << "[1/2] Partitioning into " << cfg.num_partitions
                  << " partitions  (" << p1_threads << " thread"
                  << (p1_threads > 1 ? "s" : "") << ")...\n";


    // Phase 0 (optional): pre-scan to build a load-balanced partition table.
    std::optional<RepartitionTable> rt;
    std::optional<KmcNormTable>     kmc_nt;
    std::optional<KmcPartMap>       kmc_pm;
    double t_phase0 = 0.0;

    if (cfg.strategy == PartitionStrategy::KMC) {
        if (!cfg.hide_progress)
            std::cerr << "    [0/2] Building KMC norm table (m="
                      << cfg.kmc_sig_len << ")...\n";
        kmc_nt.emplace(cfg.kmc_sig_len);

        if (!cfg.hide_progress)
            std::cerr << "    [0/2] Scanning KMC signature frequencies...\n";
        const auto t_scan = std::chrono::steady_clock::now();
        auto counts = scan_kmc_sig_counts<k>(cfg, *kmc_nt);
        kmc_pm.emplace(KmcPartMap::from_counts(*kmc_nt, counts, cfg.num_partitions));
        t_phase0 = elapsed_s(t_scan);

    } else if (cfg.strategy == PartitionStrategy::KMTRICKS) {
        if (!cfg.hide_progress)
            std::cerr << "    [0/2] Scanning minimizer frequencies (kmtricks)...\n";
        const auto t_scan = std::chrono::steady_clock::now();
        auto counts = scan_minimizer_counts<k>(cfg);
        rt.emplace(RepartitionTable::from_counts(counts, cfg.num_partitions));
        t_phase0 = elapsed_s(t_scan);
    }

    std::vector<std::ofstream> buckets(cfg.num_partitions);
    for (size_t p = 0; p < cfg.num_partitions; ++p)
        buckets[p].open(cfg.work_dir + cfg.partition_prefix() + "_" + std::to_string(p) + ".superkmers",
                        std::ios::binary);

    const auto t_part = std::chrono::steady_clock::now();
    const auto stats = cfg.strategy == PartitionStrategy::KMC
        ? partition_kmers_kmc<k, l>(cfg, buckets, *kmc_nt, *kmc_pm)
        : cfg.strategy == PartitionStrategy::KMTRICKS
            ? partition_kmers_kmtricks<k, l>(cfg, buckets, *rt)
            : partition_kmers_hash<k, l>(cfg, buckets);
    for (auto& f : buckets) f.close();

    const double t_phase1 = elapsed_s(t_part); // partitioning only (excl. pre-scan)
    if (!cfg.hide_progress)
        std::cerr << "    " << stats.seqs  << " sequences, "
                  << stats.kmers << " k-mers  ("
                  << t_phase1 << "s)\n";

    if (cfg.partition_only) {
        if (t_phase0 > 0.0)
            std::cerr << "phase0: " << t_phase0 << "s\n";
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

    const auto [total_inserted, total_written] = count_and_write<k, l>(cfg, out);

    const double t_phase2 = elapsed_s(t2);
    if (!cfg.hide_progress)
        std::cerr << "    " << total_inserted << " k-mers processed, "
                  << total_written << " unique written  ("
                  << t_phase2 << "s)\n";

    // Always emit structured phase timings for tooling (bench.sh).
    if (t_phase0 > 0.0)
        std::cerr << "phase0: " << t_phase0 << "s\n";
    std::cerr << "phase1: " << t_phase1 << "s\n"
              << "phase2: " << t_phase2 << "s\n";

    if (!cfg.hide_progress)
        std::cerr << "[done] total: " << elapsed_s(t_start) << "s\n";

    return 0;
}


// ─── Runtime dispatch (k × l) ────────────────────────────────────────────────
// kache-hash Kmer<k> encodes bases as 2 bits in a single uint64_t → k ≤ 32.
// Supported k: 9..32  (k ≥ 9 required so that at least one l in [8,k) exists)
// Supported l: 8..31  (l < k, both in [8,32])

#define TUNA_DISPATCH(K, L)  if (k == (K) && l == (L)) return run<K, L>(cfg)

inline int dispatch(uint16_t k, uint16_t l, const Config& cfg)
{
    // k = 9
    TUNA_DISPATCH( 9,  8);
    // k = 10
    TUNA_DISPATCH(10,  8);
    TUNA_DISPATCH(10,  9);
    // k = 11
    TUNA_DISPATCH(11,  8);
    TUNA_DISPATCH(11,  9);
    TUNA_DISPATCH(11, 10);
    // k = 12
    TUNA_DISPATCH(12,  8);
    TUNA_DISPATCH(12,  9);
    TUNA_DISPATCH(12, 10);
    TUNA_DISPATCH(12, 11);
    // k = 13
    TUNA_DISPATCH(13,  8);
    TUNA_DISPATCH(13,  9);
    TUNA_DISPATCH(13, 10);
    TUNA_DISPATCH(13, 11);
    TUNA_DISPATCH(13, 12);
    // k = 14
    TUNA_DISPATCH(14,  8);
    TUNA_DISPATCH(14,  9);
    TUNA_DISPATCH(14, 10);
    TUNA_DISPATCH(14, 11);
    TUNA_DISPATCH(14, 12);
    TUNA_DISPATCH(14, 13);
    // k = 15
    TUNA_DISPATCH(15,  8);
    TUNA_DISPATCH(15,  9);
    TUNA_DISPATCH(15, 10);
    TUNA_DISPATCH(15, 11);
    TUNA_DISPATCH(15, 12);
    TUNA_DISPATCH(15, 13);
    TUNA_DISPATCH(15, 14);
    // k = 16
    TUNA_DISPATCH(16,  8);
    TUNA_DISPATCH(16,  9);
    TUNA_DISPATCH(16, 10);
    TUNA_DISPATCH(16, 11);
    TUNA_DISPATCH(16, 12);
    TUNA_DISPATCH(16, 13);
    TUNA_DISPATCH(16, 14);
    TUNA_DISPATCH(16, 15);
    // k = 17
    TUNA_DISPATCH(17,  8);
    TUNA_DISPATCH(17,  9);
    TUNA_DISPATCH(17, 10);
    TUNA_DISPATCH(17, 11);
    TUNA_DISPATCH(17, 12);
    TUNA_DISPATCH(17, 13);
    TUNA_DISPATCH(17, 14);
    TUNA_DISPATCH(17, 15);
    TUNA_DISPATCH(17, 16);
    // k = 18
    TUNA_DISPATCH(18,  8);
    TUNA_DISPATCH(18,  9);
    TUNA_DISPATCH(18, 10);
    TUNA_DISPATCH(18, 11);
    TUNA_DISPATCH(18, 12);
    TUNA_DISPATCH(18, 13);
    TUNA_DISPATCH(18, 14);
    TUNA_DISPATCH(18, 15);
    TUNA_DISPATCH(18, 16);
    TUNA_DISPATCH(18, 17);
    // k = 19
    TUNA_DISPATCH(19,  8);
    TUNA_DISPATCH(19,  9);
    TUNA_DISPATCH(19, 10);
    TUNA_DISPATCH(19, 11);
    TUNA_DISPATCH(19, 12);
    TUNA_DISPATCH(19, 13);
    TUNA_DISPATCH(19, 14);
    TUNA_DISPATCH(19, 15);
    TUNA_DISPATCH(19, 16);
    TUNA_DISPATCH(19, 17);
    TUNA_DISPATCH(19, 18);
    // k = 20
    TUNA_DISPATCH(20,  8);
    TUNA_DISPATCH(20,  9);
    TUNA_DISPATCH(20, 10);
    TUNA_DISPATCH(20, 11);
    TUNA_DISPATCH(20, 12);
    TUNA_DISPATCH(20, 13);
    TUNA_DISPATCH(20, 14);
    TUNA_DISPATCH(20, 15);
    TUNA_DISPATCH(20, 16);
    TUNA_DISPATCH(20, 17);
    TUNA_DISPATCH(20, 18);
    TUNA_DISPATCH(20, 19);
    // k = 21
    TUNA_DISPATCH(21,  8);
    TUNA_DISPATCH(21,  9);
    TUNA_DISPATCH(21, 10);
    TUNA_DISPATCH(21, 11);
    TUNA_DISPATCH(21, 12);
    TUNA_DISPATCH(21, 13);
    TUNA_DISPATCH(21, 14);
    TUNA_DISPATCH(21, 15);
    TUNA_DISPATCH(21, 16);
    TUNA_DISPATCH(21, 17);
    TUNA_DISPATCH(21, 18);
    TUNA_DISPATCH(21, 19);
    TUNA_DISPATCH(21, 20);
    // k = 22
    TUNA_DISPATCH(22,  8);
    TUNA_DISPATCH(22,  9);
    TUNA_DISPATCH(22, 10);
    TUNA_DISPATCH(22, 11);
    TUNA_DISPATCH(22, 12);
    TUNA_DISPATCH(22, 13);
    TUNA_DISPATCH(22, 14);
    TUNA_DISPATCH(22, 15);
    TUNA_DISPATCH(22, 16);
    TUNA_DISPATCH(22, 17);
    TUNA_DISPATCH(22, 18);
    TUNA_DISPATCH(22, 19);
    TUNA_DISPATCH(22, 20);
    TUNA_DISPATCH(22, 21);
    // k = 23
    TUNA_DISPATCH(23,  8);
    TUNA_DISPATCH(23,  9);
    TUNA_DISPATCH(23, 10);
    TUNA_DISPATCH(23, 11);
    TUNA_DISPATCH(23, 12);
    TUNA_DISPATCH(23, 13);
    TUNA_DISPATCH(23, 14);
    TUNA_DISPATCH(23, 15);
    TUNA_DISPATCH(23, 16);
    TUNA_DISPATCH(23, 17);
    TUNA_DISPATCH(23, 18);
    TUNA_DISPATCH(23, 19);
    TUNA_DISPATCH(23, 20);
    TUNA_DISPATCH(23, 21);
    TUNA_DISPATCH(23, 22);
    // k = 24
    TUNA_DISPATCH(24,  8);
    TUNA_DISPATCH(24,  9);
    TUNA_DISPATCH(24, 10);
    TUNA_DISPATCH(24, 11);
    TUNA_DISPATCH(24, 12);
    TUNA_DISPATCH(24, 13);
    TUNA_DISPATCH(24, 14);
    TUNA_DISPATCH(24, 15);
    TUNA_DISPATCH(24, 16);
    TUNA_DISPATCH(24, 17);
    TUNA_DISPATCH(24, 18);
    TUNA_DISPATCH(24, 19);
    TUNA_DISPATCH(24, 20);
    TUNA_DISPATCH(24, 21);
    TUNA_DISPATCH(24, 22);
    TUNA_DISPATCH(24, 23);
    // k = 25
    TUNA_DISPATCH(25,  8);
    TUNA_DISPATCH(25,  9);
    TUNA_DISPATCH(25, 10);
    TUNA_DISPATCH(25, 11);
    TUNA_DISPATCH(25, 12);
    TUNA_DISPATCH(25, 13);
    TUNA_DISPATCH(25, 14);
    TUNA_DISPATCH(25, 15);
    TUNA_DISPATCH(25, 16);
    TUNA_DISPATCH(25, 17);
    TUNA_DISPATCH(25, 18);
    TUNA_DISPATCH(25, 19);
    TUNA_DISPATCH(25, 20);
    TUNA_DISPATCH(25, 21);
    TUNA_DISPATCH(25, 22);
    TUNA_DISPATCH(25, 23);
    TUNA_DISPATCH(25, 24);
    // k = 26
    TUNA_DISPATCH(26,  8);
    TUNA_DISPATCH(26,  9);
    TUNA_DISPATCH(26, 10);
    TUNA_DISPATCH(26, 11);
    TUNA_DISPATCH(26, 12);
    TUNA_DISPATCH(26, 13);
    TUNA_DISPATCH(26, 14);
    TUNA_DISPATCH(26, 15);
    TUNA_DISPATCH(26, 16);
    TUNA_DISPATCH(26, 17);
    TUNA_DISPATCH(26, 18);
    TUNA_DISPATCH(26, 19);
    TUNA_DISPATCH(26, 20);
    TUNA_DISPATCH(26, 21);
    TUNA_DISPATCH(26, 22);
    TUNA_DISPATCH(26, 23);
    TUNA_DISPATCH(26, 24);
    TUNA_DISPATCH(26, 25);
    // k = 27
    TUNA_DISPATCH(27,  8);
    TUNA_DISPATCH(27,  9);
    TUNA_DISPATCH(27, 10);
    TUNA_DISPATCH(27, 11);
    TUNA_DISPATCH(27, 12);
    TUNA_DISPATCH(27, 13);
    TUNA_DISPATCH(27, 14);
    TUNA_DISPATCH(27, 15);
    TUNA_DISPATCH(27, 16);
    TUNA_DISPATCH(27, 17);
    TUNA_DISPATCH(27, 18);
    TUNA_DISPATCH(27, 19);
    TUNA_DISPATCH(27, 20);
    TUNA_DISPATCH(27, 21);
    TUNA_DISPATCH(27, 22);
    TUNA_DISPATCH(27, 23);
    TUNA_DISPATCH(27, 24);
    TUNA_DISPATCH(27, 25);
    TUNA_DISPATCH(27, 26);
    // k = 28
    TUNA_DISPATCH(28,  8);
    TUNA_DISPATCH(28,  9);
    TUNA_DISPATCH(28, 10);
    TUNA_DISPATCH(28, 11);
    TUNA_DISPATCH(28, 12);
    TUNA_DISPATCH(28, 13);
    TUNA_DISPATCH(28, 14);
    TUNA_DISPATCH(28, 15);
    TUNA_DISPATCH(28, 16);
    TUNA_DISPATCH(28, 17);
    TUNA_DISPATCH(28, 18);
    TUNA_DISPATCH(28, 19);
    TUNA_DISPATCH(28, 20);
    TUNA_DISPATCH(28, 21);
    TUNA_DISPATCH(28, 22);
    TUNA_DISPATCH(28, 23);
    TUNA_DISPATCH(28, 24);
    TUNA_DISPATCH(28, 25);
    TUNA_DISPATCH(28, 26);
    TUNA_DISPATCH(28, 27);
    // k = 29
    TUNA_DISPATCH(29,  8);
    TUNA_DISPATCH(29,  9);
    TUNA_DISPATCH(29, 10);
    TUNA_DISPATCH(29, 11);
    TUNA_DISPATCH(29, 12);
    TUNA_DISPATCH(29, 13);
    TUNA_DISPATCH(29, 14);
    TUNA_DISPATCH(29, 15);
    TUNA_DISPATCH(29, 16);
    TUNA_DISPATCH(29, 17);
    TUNA_DISPATCH(29, 18);
    TUNA_DISPATCH(29, 19);
    TUNA_DISPATCH(29, 20);
    TUNA_DISPATCH(29, 21);
    TUNA_DISPATCH(29, 22);
    TUNA_DISPATCH(29, 23);
    TUNA_DISPATCH(29, 24);
    TUNA_DISPATCH(29, 25);
    TUNA_DISPATCH(29, 26);
    TUNA_DISPATCH(29, 27);
    TUNA_DISPATCH(29, 28);
    // k = 30
    TUNA_DISPATCH(30,  8);
    TUNA_DISPATCH(30,  9);
    TUNA_DISPATCH(30, 10);
    TUNA_DISPATCH(30, 11);
    TUNA_DISPATCH(30, 12);
    TUNA_DISPATCH(30, 13);
    TUNA_DISPATCH(30, 14);
    TUNA_DISPATCH(30, 15);
    TUNA_DISPATCH(30, 16);
    TUNA_DISPATCH(30, 17);
    TUNA_DISPATCH(30, 18);
    TUNA_DISPATCH(30, 19);
    TUNA_DISPATCH(30, 20);
    TUNA_DISPATCH(30, 21);
    TUNA_DISPATCH(30, 22);
    TUNA_DISPATCH(30, 23);
    TUNA_DISPATCH(30, 24);
    TUNA_DISPATCH(30, 25);
    TUNA_DISPATCH(30, 26);
    TUNA_DISPATCH(30, 27);
    TUNA_DISPATCH(30, 28);
    TUNA_DISPATCH(30, 29);
    // k = 31
    TUNA_DISPATCH(31,  8);
    TUNA_DISPATCH(31,  9);
    TUNA_DISPATCH(31, 10);
    TUNA_DISPATCH(31, 11);
    TUNA_DISPATCH(31, 12);
    TUNA_DISPATCH(31, 13);
    TUNA_DISPATCH(31, 14);
    TUNA_DISPATCH(31, 15);
    TUNA_DISPATCH(31, 16);
    TUNA_DISPATCH(31, 17);
    TUNA_DISPATCH(31, 18);
    TUNA_DISPATCH(31, 19);
    TUNA_DISPATCH(31, 20);
    TUNA_DISPATCH(31, 21);
    TUNA_DISPATCH(31, 22);
    TUNA_DISPATCH(31, 23);
    TUNA_DISPATCH(31, 24);
    TUNA_DISPATCH(31, 25);
    TUNA_DISPATCH(31, 26);
    TUNA_DISPATCH(31, 27);
    TUNA_DISPATCH(31, 28);
    TUNA_DISPATCH(31, 29);
    TUNA_DISPATCH(31, 30);
    // k = 32
    TUNA_DISPATCH(32,  8);
    TUNA_DISPATCH(32,  9);
    TUNA_DISPATCH(32, 10);
    TUNA_DISPATCH(32, 11);
    TUNA_DISPATCH(32, 12);
    TUNA_DISPATCH(32, 13);
    TUNA_DISPATCH(32, 14);
    TUNA_DISPATCH(32, 15);
    TUNA_DISPATCH(32, 16);
    TUNA_DISPATCH(32, 17);
    TUNA_DISPATCH(32, 18);
    TUNA_DISPATCH(32, 19);
    TUNA_DISPATCH(32, 20);
    TUNA_DISPATCH(32, 21);
    TUNA_DISPATCH(32, 22);
    TUNA_DISPATCH(32, 23);
    TUNA_DISPATCH(32, 24);
    TUNA_DISPATCH(32, 25);
    TUNA_DISPATCH(32, 26);
    TUNA_DISPATCH(32, 27);
    TUNA_DISPATCH(32, 28);
    TUNA_DISPATCH(32, 29);
    TUNA_DISPATCH(32, 30);
    TUNA_DISPATCH(32, 31);

    std::cerr << "tuna: error: unsupported combination k=" << k << " l=" << l << "\n"
              << "  k must be in [9,32]  (kache-hash Kmer fits in 64 bits, k ≤ 32)\n"
              << "  l must be in [8,k-1] and strictly less than k\n";
    return 1;
}

#undef TUNA_DISPATCH
