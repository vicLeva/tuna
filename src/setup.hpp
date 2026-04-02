#pragma once

// Partition-count auto-tuning utilities.
// Shared between the CLI (main.cpp) and the public library API (include/tuna/tuna.hpp).

#include "Config.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <sys/resource.h>

// Read a cache level's size from Linux sysfs (cpu0, shared or per-core).
// Scans index entries for the requested level. Falls back to `fallback` if unavailable.
inline size_t sysfs_cache_size(int target_level, size_t fallback)
{
    for (int idx = 0; idx <= 5; ++idx) {
        const std::string base =
            "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(idx);
        {
            std::ifstream lf(base + "/level");
            if (!lf) continue;
            int level = 0;
            if (!(lf >> level) || level != target_level) continue;
        }
        // For L3 we want "Unified"; skip instruction-only entries.
        {
            std::ifstream tf(base + "/type");
            if (tf) {
                std::string type;
                if (tf >> type && type == "Instruction") continue;
            }
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
    return fallback;
}

inline size_t l2_cache_per_core() { return sysfs_cache_size(2, 256u << 10); }
inline size_t l3_cache_total()    { return sysfs_cache_size(3, 30u  << 20); }


// Round n up to the next power of 2 (n = 0 or 1 → 1).
inline uint32_t round_pow2(uint32_t n)
{
    if (n <= 1) return 1;
    if ((n & (n - 1)) == 0) return n;
    while (n & (n - 1)) n &= n - 1;
    return n << 1;
}


// Auto-tune partition count for phase-2 LLC residency.
//
// Goal: all t active phase-2 tables fit in L3 cache simultaneously, so that
// k-mer upserts hit L3 rather than DRAM.
//
//   n_target = ceil_pow2( t × est_unique × bytes_per_slot / (L3_USABLE) )
//
// where:
//   est_unique       — estimated unique k-mers (derived from input bytes)
//   bytes_per_slot   — flat_t + Metadata per slot, adjusted for k and load factor
//   L3_USABLE        — 75% of L3 (leave headroom for OS + phase-1 buffers)
//
// For multi-file runs where k-mers overlap across files (e.g. 200 E. coli
// genomes), total bytes greatly overstates unique k-mers.  We cap the estimate
// at the largest single-file contribution to avoid over-partitioning.
//
// writer_header_bytes — sizeof(SuperkmerWriter); passed by caller to keep this
//   header free of the heavy partition_hash.hpp include.
// k                  — k-mer length; determines sizeof(flat_t).
// num_threads        — number of phase-2 worker threads (= active tables).
inline uint32_t auto_tune_partitions(
    const std::vector<std::string>& input_files,
    size_t                          writer_header_bytes,
    uint16_t                        k,
    uint32_t                        num_threads)
{
    namespace fs = std::filesystem;

    // ── FD budget ──────────────────────────────────────────────────────────
    struct rlimit rl{};
    size_t fd_limit = 1024;
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
        fd_limit = static_cast<size_t>(rl.rlim_cur);
    const size_t max_parts = (fd_limit > 32) ? fd_limit - 32 : 16;

    // ── Writer struct cap (phase 1) ────────────────────────────────────────
    // Largest power of 2 s.t. n × writer_header_bytes ≤ L2 per core.
    const size_t l2    = l2_cache_per_core();
    const size_t ratio = l2 / writer_header_bytes;
    size_t writer_cap  = 4096;
    while ((writer_cap << 1) <= ratio) writer_cap <<= 1;
    const size_t WRITER_CACHE_CAP = std::min(writer_cap, size_t(32768));

    // ── Estimate uncompressed input ────────────────────────────────────────
    constexpr uint64_t GZ_EXPAND = 6;
    uint64_t total_effective = 0;
    uint64_t max_single      = 0;   // largest single file's contribution
    for (const auto& f : input_files) {
        std::error_code ec;
        uint64_t sz = fs::file_size(f, ec);
        if (ec) continue;
        if (f.size() >= 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
            sz *= GZ_EXPAND;
        total_effective += sz;
        max_single = std::max(max_single, sz);
    }

    // ── Estimate unique k-mers ─────────────────────────────────────────────
    // For multi-file same-species runs (avg file < 50 MB), k-mers overlap
    // heavily — unique count is closer to the single largest file than the sum.
    // For large or few files, use the total (each file is likely a different sample).
    const size_t n_files = input_files.size();
    uint64_t est_unique_bytes;
    if (n_files > 1) {
        const uint64_t avg = total_effective / n_files;
        constexpr uint64_t SMALL_FILE_THRESHOLD = 50ULL << 20;   // 50 MB
        est_unique_bytes = (avg < SMALL_FILE_THRESHOLD)
            ? max_single           // e.g. 200 E. coli: use one genome's worth
            : total_effective;     // e.g. 10 human samples: all are unique
    } else {
        est_unique_bytes = total_effective;
    }

    // Estimate unique k-mers from uncompressed FASTA bytes.
    //
    // In a genomic FASTA file ~90% of bytes are sequence characters, of which
    // ~85% are ACTG (the rest are N or line endings).  That gives ~76% of file
    // bytes as countable bases.
    //
    // For genomic data, the number of *unique* k-mers is approximately equal to
    // the number of bases: nearly every position yields a distinct k-mer
    // (especially for k≥21).  Dividing by k is wrong — it gives total k-mers
    // but unique ≈ total for typical genomes.
    const uint64_t est_unique_kmers = est_unique_bytes * 76 / 100;

    // ── Per-slot memory in the hash table ─────────────────────────────────
    // flat_t = Kmer<k> + uint32_t count.
    // Kmer<k>: ceil(k/32) × 8 bytes (one uint64_t per 32 bases).
    // Metadata: cs (1 B) + min_coord (1 B) per slot = 2 B/slot.
    // At load factor 0.8: each unique k-mer occupies 1/0.8 = 1.25 slots.
    const size_t kmer_words     = (k + 31) / 32;
    const size_t flat_t_bytes   = kmer_words * 8 + 4;   // Kmer<k> + uint32_t
    const size_t meta_bytes     = 2;                     // cs + min_coord
    const size_t bytes_per_slot = flat_t_bytes + meta_bytes;
    // bytes per unique k-mer (at lf=0.8, 1.25 slots/unique):
    const size_t bytes_per_kmer = (bytes_per_slot * 5 + 3) / 4;  // × 1.25, integer

    // ── L3-residency target ────────────────────────────────────────────────
    // We want: num_threads × table_bytes < 0.75 × L3
    // table_bytes = est_unique_kmers_per_part × bytes_per_kmer
    //            = (est_unique_kmers / n) × bytes_per_kmer
    // Solving for n:
    //   n > num_threads × est_unique_kmers × bytes_per_kmer / (0.75 × L3)
    const size_t l3        = l3_cache_total();
    const size_t l3_usable = l3 * 3 / 4;
    const size_t t         = std::max(uint32_t(1), num_threads);

    // Compute n_l3 = ceil(numerator / l3_usable), then round up to pow2.
    const uint64_t numerator = static_cast<uint64_t>(t)
                             * est_unique_kmers
                             * bytes_per_kmer;
    const size_t n_l3_raw = static_cast<size_t>(
        (numerator + l3_usable - 1) / l3_usable);
    size_t n_l3 = 16;
    while (n_l3 < n_l3_raw) n_l3 <<= 1;

    // ── Throughput floor ───────────────────────────────────────────────────
    // Even when tables fit in L3 at low n, more partitions = smaller work
    // units per round = better thread utilisation and fewer resizes.
    // Use total_effective / 2 MB as a throughput-oriented lower bound
    // (the original heuristic), uncapped by the multi-file correction above.
    const size_t n_tp_raw = static_cast<size_t>(
        std::max(uint64_t(1), total_effective >> 21));   // / 2 MB
    size_t n_throughput = 16;
    while (n_throughput < n_tp_raw) n_throughput <<= 1;

    // Final n: whichever constraint is more demanding.
    size_t n = std::max(n_l3, n_throughput);

    // ── Apply caps ─────────────────────────────────────────────────────────
    n = std::min(n, std::min(max_parts, WRITER_CACHE_CAP));
    n = std::max(n, size_t(16));

    return static_cast<uint32_t>(n);
}
