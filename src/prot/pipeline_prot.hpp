#pragma once

// Protein pipeline orchestrator.
// Selects in-memory vs disk mode based on estimated superkmer size vs available RAM.
// Emits structured timing on stderr: "phase1: Xs\nphase2: Xs\n"

#include "config_prot.hpp"
#include "partition_prot.hpp"
#include "count_prot.hpp"
#include "setup.hpp"   // auto_tune_partitions, l3_cache_total

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef __linux__
#  include <sys/sysinfo.h>
#endif

namespace prot {

constexpr uint64_t GZ_EXPAND_PROT = 6;

inline uint64_t available_ram_prot() {
#ifdef __linux__
    struct sysinfo si{};
    if (sysinfo(&si) == 0)
        return static_cast<uint64_t>(si.freeram + si.bufferram) * si.mem_unit;
#endif
    return 0;
}

inline double elapsed_s_prot(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

inline std::string fmt_s_prot(double s) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.3fs", s);
    return buf;
}


template <uint16_t k, uint16_t m>
int run_prot(const ProtConfig& cfg) {
    // ── Resolve output stream ──────────────────────────────────────────────
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!cfg.output_path.empty() && cfg.output_path != "-") {
        out_file.open(cfg.output_path);
        if (!out_file) throw std::runtime_error("Cannot open output: " + cfg.output_path);
        out = &out_file;
    }

    // ── Auto-tune partition count ──────────────────────────────────────────
    uint32_t n_parts = cfg.num_partitions;
    if (n_parts == 0) {
        // For proteins, use a simplified heuristic: total input bytes / 1 MB,
        // clamped to [16, 4096], rounded to power of 2.
        uint64_t total = 0;
        for (const auto& f : cfg.input_files) {
            std::error_code ec;
            uint64_t fsz = std::filesystem::file_size(f, ec);
            if (ec) continue;
            const bool gz = f.size() > 3 && f.compare(f.size()-3,3,".gz")==0;
            total += gz ? fsz * GZ_EXPAND_PROT : fsz;
        }
        n_parts = 16;
        const uint64_t target = total >> 20;  // / 1 MB
        while (n_parts < target && n_parts < 4096) n_parts <<= 1;
        n_parts = std::max(n_parts, uint32_t(16));
    }
    // Ensure power of 2.
    {
        uint32_t p = 1;
        while (p < n_parts) p <<= 1;
        n_parts = p;
    }

    // ── RAM budget and pipeline mode ───────────────────────────────────────
    uint64_t est_packed = 0;
    for (const auto& f : cfg.input_files) {
        std::error_code ec;
        uint64_t fsz = std::filesystem::file_size(f, ec);
        if (ec) continue;
        const bool gz = f.size() > 3 && f.compare(f.size()-3,3,".gz")==0;
        // Protein seqs: ~5/8 packing ratio of uncompressed bytes.
        est_packed += gz
            ? static_cast<uint64_t>(fsz * GZ_EXPAND_PROT * 5 / 8)
            : fsz * 5 / 8;
    }

    const uint64_t avail = cfg.ram_budget_bytes > 0
        ? cfg.ram_budget_bytes : available_ram_prot();
    const bool use_mem = avail > 0 && est_packed < avail * 6 / 10;

    const size_t disk_write_budget = [&]() -> size_t {
        if (avail == 0) return size_t(64) << 20;
        const uint64_t bgt = avail * 3 / 10 / cfg.num_threads;
        return static_cast<size_t>(std::min(bgt, uint64_t(512) << 20));
    }();

    if (!cfg.hide_progress) {
        std::cerr << "tuna-prot k=" << cfg.k << " m=" << cfg.m
                  << " n=" << n_parts << " t=" << cfg.num_threads
                  << (use_mem ? " [in-memory]" : " [disk]") << "\n";
    }

    // ── Phase 1 ────────────────────────────────────────────────────────────
    const auto t1 = std::chrono::steady_clock::now();

    if (use_mem) {
        std::vector<std::string> part_bufs(n_parts);
        [[maybe_unused]] const PartitionStats ps = partition_kmers_mem<k, m>(
            cfg, part_bufs);

        std::cerr << "phase1: " << fmt_s_prot(elapsed_s_prot(t1)) << "\n";

        if (cfg.stop_after_phase1) return 0;

        // ── Phase 2+3 ──────────────────────────────────────────────────────
        const auto t2 = std::chrono::steady_clock::now();
        count_and_write_prot_mem<k, m>(cfg, part_bufs, *out);
        std::cerr << "phase2: " << fmt_s_prot(elapsed_s_prot(t2)) << "\n";

    } else {
        // Disk mode: open partition files.
        std::vector<std::ofstream> buckets(n_parts);
        for (size_t p = 0; p < n_parts; ++p) {
            const std::string path = prot_partition_path(cfg.work_dir, p);
            buckets[p].open(path, std::ios::binary | std::ios::trunc);
            if (!buckets[p]) throw std::runtime_error("Cannot create: " + path);
        }

        [[maybe_unused]] const PartitionStats ps = partition_kmers_disk<k, m>(
            cfg, cfg.work_dir, buckets, disk_write_budget);
        for (auto& b : buckets) b.close();

        std::cerr << "phase1: " << fmt_s_prot(elapsed_s_prot(t1)) << "\n";

        if (cfg.stop_after_phase1) return 0;

        const auto t2 = std::chrono::steady_clock::now();
        count_and_write_prot_disk<k, m>(cfg, cfg.work_dir, n_parts, *out);
        std::cerr << "phase2: " << fmt_s_prot(elapsed_s_prot(t2)) << "\n";

        // Clean up partition files.
        if (!cfg.keep_tmp) {
            for (size_t p = 0; p < n_parts; ++p) {
                std::error_code ec;
                std::filesystem::remove(prot_partition_path(cfg.work_dir, p), ec);
            }
        }
    }

    return 0;
}


// ── k/m dispatch ─────────────────────────────────────────────────────────────
// Template instantiation for k in [4,12] and m in [2, k-1].
// Uses a nested switch to select the correct (k, m) pair at runtime.

inline int dispatch_run_prot(const ProtConfig& cfg) {
#define PROT_DISPATCH_M(K)                                              \
    switch (cfg.m) {                                                    \
        case  2: return run_prot<K,  2>(cfg);                          \
        case  3: return run_prot<K,  3>(cfg);                          \
        case  4: if constexpr(K >  4) return run_prot<K,  4>(cfg); break; \
        case  5: if constexpr(K >  5) return run_prot<K,  5>(cfg); break; \
        case  6: if constexpr(K >  6) return run_prot<K,  6>(cfg); break; \
        case  7: if constexpr(K >  7) return run_prot<K,  7>(cfg); break; \
        case  8: if constexpr(K >  8) return run_prot<K,  8>(cfg); break; \
        case  9: if constexpr(K >  9) return run_prot<K,  9>(cfg); break; \
        case 10: if constexpr(K > 10) return run_prot<K, 10>(cfg); break; \
        case 11: if constexpr(K > 11) return run_prot<K, 11>(cfg); break; \
        default: break;                                                 \
    }

    switch (cfg.k) {
        case  4: PROT_DISPATCH_M( 4) break;
        case  5: PROT_DISPATCH_M( 5) break;
        case  6: PROT_DISPATCH_M( 6) break;
        case  7: PROT_DISPATCH_M( 7) break;
        case  8: PROT_DISPATCH_M( 8) break;
        case  9: PROT_DISPATCH_M( 9) break;
        case 10: PROT_DISPATCH_M(10) break;
        case 11: PROT_DISPATCH_M(11) break;
        case 12: PROT_DISPATCH_M(12) break;
        default: break;
    }
#undef PROT_DISPATCH_M
    throw std::invalid_argument(
        "Unsupported k=" + std::to_string(cfg.k) +
        " m=" + std::to_string(cfg.m) +
        " (k in [4,12], m in [2,k-1])");
}

} // namespace prot
