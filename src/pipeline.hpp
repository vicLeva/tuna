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

    const auto t1 = std::chrono::steady_clock::now();

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
// Supported k: 17 21 27 31
// Supported l: 11 13 15 17  (with l < k)

#define TUNA_DISPATCH(K, L)  if (k == (K) && l == (L)) return run<K, L>(cfg)

inline int dispatch(uint16_t k, uint16_t l, const Config& cfg)
{
    // k = 17
    TUNA_DISPATCH(17, 11);
    TUNA_DISPATCH(17, 13);
    TUNA_DISPATCH(17, 15);
    // k = 21
    TUNA_DISPATCH(21, 11);
    TUNA_DISPATCH(21, 13);
    TUNA_DISPATCH(21, 15);
    TUNA_DISPATCH(21, 17);
    // k = 27
    TUNA_DISPATCH(27, 11);
    TUNA_DISPATCH(27, 13);
    TUNA_DISPATCH(27, 15);
    TUNA_DISPATCH(27, 17);
    // k = 31
    TUNA_DISPATCH(31, 11);
    TUNA_DISPATCH(31, 13);
    TUNA_DISPATCH(31, 15);
    TUNA_DISPATCH(31, 17);   // default

    std::cerr << "tuna: error: unsupported combination k=" << k << " l=" << l << "\n"
              << "  k must be in { 17, 21, 27, 31 }  (kache-hash Kmer fits in 64 bits, k ≤ 32)\n"
              << "  l must be in { 11, 13, 15, 17 } and strictly less than k\n";
    return 1;
}

#undef TUNA_DISPATCH
