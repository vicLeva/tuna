#pragma once

// Partition-count auto-tuning utilities.
// Shared between the CLI (main.cpp) and the public library API (include/tuna/tuna.hpp).

#include "Config.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <sys/resource.h>

// Read L2 data cache size per core from Linux sysfs.
// Scans cache index entries for level=2 to avoid assuming a fixed index number.
// Falls back to 256 KB if unavailable.
inline size_t l2_cache_per_core()
{
    for (int idx = 0; idx <= 4; ++idx) {
        const std::string base =
            "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(idx);
        {
            std::ifstream lf(base + "/level");
            if (!lf) continue;
            int level = 0;
            if (!(lf >> level) || level != 2) continue;
        }
        std::ifstream sf(base + "/size");
        if (!sf) continue;
        std::string s;
        if (!(sf >> s) || s.empty()) continue;
        char* end;
        size_t val = std::strtoul(s.c_str(), &end, 10);
        if      (*end == 'K' || *end == 'k') val <<= 10;
        else if (*end == 'M' || *end == 'm') val <<= 20;
        if (val) return val;
    }
    return 256u << 10;
}


// Round n up to the next power of 2 (n = 0 or 1 → 1).
inline uint32_t round_pow2(uint32_t n)
{
    if (n <= 1) return 1;
    if ((n & (n - 1)) == 0) return n;
    while (n & (n - 1)) n &= n - 1;
    return n << 1;
}


// Auto-tune partition count to ~2 MB estimated uncompressed sequence per partition.
//
// writer_header_bytes — sizeof(SuperkmerWriter); passed by the caller to keep this
//   header lightweight (no need to include the heavy partition_hash.hpp here).
//
// The writer cache cap ensures Phase 1's writer array fits in per-core L2 so the
// minimizer-dispatch loop does not thrash the cache.
inline uint32_t auto_tune_partitions(
    const std::vector<std::string>& input_files,
    size_t                          writer_header_bytes)
{
    namespace fs = std::filesystem;

    // ── FD budget ──────────────────────────────────────────────────────────
    struct rlimit rl{};
    size_t fd_limit = 1024;
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
        fd_limit = static_cast<size_t>(rl.rlim_cur);
    const size_t max_parts = (fd_limit > 32) ? fd_limit - 32 : 16;

    // ── Writer cache cap ───────────────────────────────────────────────────
    // Largest power of 2 s.t. n × writer_header_bytes ≤ L2 per core.
    const size_t l2    = l2_cache_per_core();
    const size_t ratio = l2 / writer_header_bytes;
    size_t writer_cap  = 4096;
    while ((writer_cap << 1) <= ratio) writer_cap <<= 1;
    const size_t WRITER_CACHE_CAP = std::min(writer_cap, size_t(32768));

    // ── Estimate total uncompressed input ──────────────────────────────────
    constexpr uint64_t GZ_EXPAND = 6;   // typical expansion factor for gz genomic data
    uint64_t total_effective = 0;
    for (const auto& f : input_files) {
        std::error_code ec;
        uint64_t sz = fs::file_size(f, ec);
        if (ec) continue;
        if (f.size() >= 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
            sz *= GZ_EXPAND;
        total_effective += sz;
    }

    // ── Base n: next_pow2(total / 2 MB), clamped ──────────────────────────
    const size_t raw = static_cast<size_t>(
        std::max(uint64_t(1), total_effective >> 21));   // / 2 MB
    size_t n = 1;
    while (n < raw) n <<= 1;
    n = std::max(n, size_t(16));
    n = std::min(n, std::min(max_parts, WRITER_CACHE_CAP));

    // ── Small-file cap ─────────────────────────────────────────────────────
    // For many small equal-size files (e.g. 200 E. coli genomes), cap at
    // next_pow2(n_files) to avoid unnecessary partitions.
    const size_t n_files = input_files.size();
    if (n_files > 1) {
        const uint64_t avg_effective = total_effective / n_files;
        constexpr uint64_t SMALL_FILE_THRESHOLD = 50ULL << 20;   // 50 MB
        if (avg_effective < SMALL_FILE_THRESHOLD) {
            size_t file_cap = 1;
            while (file_cap < n_files) file_cap <<= 1;
            n = std::min(n, file_cap);
        }
    }

    return static_cast<uint32_t>(n);
}
