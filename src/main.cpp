// tuna — streaming k-mer counter built on helicase + kache-hash
//
// Pipeline:
//   Phase 1 — stream FASTA/Q through a minimizer iterator, cut superkmers at
//             minimizer-partition boundaries, write binary superkmer records
//             to per-partition files on disk.
//   Phase 2 — replay each partition through a Kmer_Window, upsert every k-mer
//             into a Streaming_Kmer_Hash_Table with increment semantics.
//   Output  — iterate the table, apply ci/cx filters, write <kmer>\t<count>.

#include "CLI.hpp"
#include "pipeline.hpp"
#include "superkmer_io.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sys/resource.h>

// Read L2 cache size per core from Linux sysfs.
// Scans index0..index4 for level=2 to avoid assuming a fixed index number —
// on some CPUs (AMD EPYC, unified-L1 Intel) index2 is L3, not L2.
// Falls back to 256 KB if no level-2 entry is found.
static size_t l2_cache_per_core()
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


// Auto-tune the partition count when -n is not given (cfg.num_partitions == 0).
//
// Target: ~2 MB of estimated uncompressed sequence per partition.
// GZ files are expanded by GZ_EXPAND to approximate uncompressed size.
//
// Two caps are applied (in addition to the fd-limit budget):
//
//   WRITER_CACHE_CAP — keeps the writer-header array (n×sizeof(SuperkmerWriter))
//     within the per-core L2 cache.  Phase 1 performance ∝ n_parts when the
//     array exceeds L2 (cache thrashing on random-access appends):
//       tara 5.9 GB gz @ n=32768 → phase1 ≈ 387 s
//       tara 5.9 GB gz @ n= 8192 → phase1 ≈  89 s  (4× fewer → 4× faster)
//     Cap = prev_pow2(L2_bytes / sizeof(SuperkmerWriter)), clamped to [4096, 32768].
//     Detected at runtime from sysfs; falls back to 4096 (256 KB L2).
//
//   small-file cap — for multi-file runs where average file < 50 MB: cap n at
//     next_pow2(n_files) to avoid over-partitioning redundant genomes.
//
// n is rounded to the next power of 2 and clamped to [16, min(caps)].
static uint32_t auto_tune_partitions(const std::vector<std::string>& input_files)
{
    namespace fs = std::filesystem;

    // ── FD budget ─────────────────────────────────────────────────────────
    struct rlimit rl{};
    size_t fd_limit = 1024;
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
        fd_limit = static_cast<size_t>(rl.rlim_cur);
    const size_t max_parts = (fd_limit > 32) ? fd_limit - 32 : 16;

    // ── Writer cache cap ──────────────────────────────────────────────────
    // Largest power of 2 s.t. n×sizeof(SuperkmerWriter) ≤ L2 per core.
    //   L2=256 KB → cap=4096  (server core, validated: tara phase1=89s)
    //   L2=512 KB → cap=8192
    //   L2≥2  MB  → cap=32768 (clamped)
    constexpr size_t WRITER_HEADER_BYTES = sizeof(SuperkmerWriter);
    static_assert(WRITER_HEADER_BYTES >= 8 && WRITER_HEADER_BYTES <= 64,
                  "SuperkmerWriter size out of expected range");
    const size_t l2    = l2_cache_per_core();
    const size_t ratio = l2 / WRITER_HEADER_BYTES;
    size_t writer_cap  = 4096;
    while ((writer_cap << 1) <= ratio) writer_cap <<= 1;
    const size_t WRITER_CACHE_CAP = std::min(writer_cap, size_t(32768));

    // ── Estimate total uncompressed input ─────────────────────────────────
    uint64_t total_effective = 0;
    for (const auto& f : input_files) {
        std::error_code ec;
        uint64_t sz = fs::file_size(f, ec);
        if (f.size() >= 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
            sz *= GZ_EXPAND;
        total_effective += sz;
    }

    // ── Base n: next_pow2(total / 2MB), clamped ───────────────────────────
    const size_t raw = static_cast<size_t>(
        std::max(uint64_t(1), total_effective >> 21));  // / 2 MB
    size_t n = 1;
    while (n < raw) n <<= 1;
    n = std::max(n, size_t(16));
    n = std::min(n, std::min(max_parts, WRITER_CACHE_CAP));

    // ── Small-file cap ────────────────────────────────────────────────────
    // For many small files, total_effective overestimates unique k-mers
    // (redundant genomes): cap at next_pow2(n_files) to avoid excess partitions.
    const size_t n_files = input_files.size();
    if (n_files > 1) {
        const uint64_t avg_effective = total_effective / n_files;
        constexpr uint64_t SMALL_FILE_THRESHOLD = 50ULL << 20;  // 50 MB
        if (avg_effective < SMALL_FILE_THRESHOLD) {
            size_t file_cap = 1;
            while (file_cap < n_files) file_cap <<= 1;
            n = std::min(n, file_cap);
        }
    }

    return static_cast<uint32_t>(n);
}


int main(int argc, char* argv[])
{
    Config cfg;
    if (!parse_args(argc, argv, cfg)) {
        print_usage(argv[0]);
        return 1;
    }

    namespace fs = std::filesystem;

    // ── Raise soft FD limit to hard limit (best-effort) ───────────────────
    {
        struct rlimit rl{};
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0) {
            rl.rlim_cur = rl.rlim_max;
            setrlimit(RLIMIT_NOFILE, &rl);
        }
    }

    if (cfg.num_partitions == 0)
        cfg.num_partitions = auto_tune_partitions(cfg.input_files);

    // Round up to next power of 2 so partition_fn can use bitmask instead of division.
    {
        uint32_t n = cfg.num_partitions;
        if (n > 1 && (n & (n - 1)) != 0) {
            while (n & (n - 1)) n &= n - 1;
            n <<= 1;
            cfg.num_partitions = n;
        }
    }

    // ── Working directory ──────────────────────────────────────────────────
    bool own_work_dir = false;
    if (cfg.work_dir.empty()) {
        fs::path out_parent = fs::path(cfg.output_file).parent_path();
        if (out_parent.empty()) out_parent = fs::current_path();
        cfg.work_dir = (out_parent / "tuna_tmp").string();
        own_work_dir = true;
    }
    if (cfg.work_dir.back() != '/') cfg.work_dir += '/';
    fs::create_directories(cfg.work_dir);

    if (!cfg.hide_progress) {
        std::cerr << "tuna  k=" << cfg.k
                  << "  m=" << cfg.l
                  << "  n=" << cfg.num_partitions
                  << "  t=" << cfg.num_threads;
        if (cfg.ci > 1 || cfg.cx != std::numeric_limits<uint64_t>::max())
            std::cerr << "  ci=" << cfg.ci << "  cx=" << cfg.cx;
        std::cerr << "\n";
    }

    const int rc = dispatch(cfg.k, cfg.l, cfg);

    // ── Cleanup ────────────────────────────────────────────────────────────
    if (!cfg.keep_tmp) {
        if (own_work_dir)
            fs::remove_all(cfg.work_dir);
        else {
            // Only remove the partition files we created, not the whole directory.
            for (size_t p = 0; p < cfg.num_partitions; ++p)
                fs::remove(partition_path(cfg.work_dir, p));
        }
    }

    return rc;
}
