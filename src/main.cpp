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

    // ── Auto-tune partition count from total input size ────────────────────
    // Triggered when -n is not given (num_partitions == 0).
    // Target: ~2 MB of estimated uncompressed sequence per partition.
    //
    // .gz files are expanded by GZ_EXPAND to approximate uncompressed size.
    //
    // Phase 1 performance is bounded by n_parts: each consumer thread keeps one
    // SuperkmerWriter per partition, so writer-buffer cache footprint = n_parts×4 KB.
    // Benchmarks show phase1 ∝ n_parts (writer buffers exceed L3 → cache thrashing):
    //   tara 5.9 GB gz @ n=32768 → phase1 ≈ 387 s
    //   tara 5.9 GB gz @ n=16384 → phase1 ≈ 202 s  (2× fewer → 2× faster)
    // GZ_EXPAND=6 pushed tara to n=32768; we cap at 8192 so that writer buffers
    // (8192×4 KB = 32 MB) stay within a typical server L3 cache.  This cap also
    // corresponds to writer headers (8192×24 B = 192 KB) fitting in the L2 cache,
    // which is the regime where random-access append is fast.
    // Phase 2 tables grow proportionally but the cross-superkmer prefetch hides
    // the resulting L3 misses (empirically: <10 % phase2 slowdown per 2× n_parts
    // reduction) so the net effect is a large improvement.
    //
    // n is rounded to the next power of 2 and clamped to [16, min(WRITER_CAP, fd_budget)].
    if (cfg.num_partitions == 0) {
        struct rlimit rl{};
        size_t fd_limit = 1024;
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
            fd_limit = static_cast<size_t>(rl.rlim_cur);
        const size_t max_parts = (fd_limit > 32) ? fd_limit - 32 : 16;

        // WRITER_CACHE_CAP: largest power of 2 such that the writer header array
        // (n_parts × sizeof(SuperkmerWriter) ≈ 40 B) just fits within the per-core
        // L2 cache, keeping random-access appends to writers[p] out of L3/RAM.
        //   L2=256 KB → 6554 → next_pow2 = 8192  (typical server core)
        //   L2=512 KB → 13107 → next_pow2 = 16384 (workstation core)
        //   L2≥2 MB   → capped at 32768
        // Detected at runtime from /sys/devices/system/cpu/cpu0/cache/index2/size;
        // falls back to 8192 (256 KB L2 assumption) if unavailable.
        // sizeof(SuperkmerWriter) = sizeof(std::string) + sizeof(size_t).
        // On GCC/libstdc++ x86-64: 32 + 8 = 40 B.  The static_assert below
        // catches mismatches on other toolchains at compile time.
        constexpr size_t WRITER_HEADER_BYTES = sizeof(SuperkmerWriter);
        static_assert(WRITER_HEADER_BYTES >= 8 && WRITER_HEADER_BYTES <= 64,
                      "SuperkmerWriter size out of expected range — update WRITER_HEADER_BYTES");
        const size_t l2     = l2_cache_per_core();
        const size_t ratio  = l2 / WRITER_HEADER_BYTES;
        // Largest power of 2 such that n × WRITER_HEADER_BYTES ≤ L2 (prev_pow2).
        // Using next_pow2 would give a cap 2× larger than the L2 budget.
        size_t writer_cache_cap = 4096;
        while ((writer_cache_cap << 1) <= ratio) writer_cache_cap <<= 1;
        const size_t WRITER_CACHE_CAP = std::min(writer_cache_cap, size_t(32768));

        constexpr uint64_t GZ_EXPAND = 6;  // typical genomic FASTA/FASTQ compression ratio

        uint64_t total_effective = 0;
        for (const auto& f : cfg.input_files) {
            std::error_code ec;
            uint64_t sz = fs::file_size(f, ec);
            if (f.size() >= 3 && f.compare(f.size() - 3, 3, ".gz") == 0)
                sz *= GZ_EXPAND;
            total_effective += sz;
        }
        const size_t raw = static_cast<size_t>(
            std::max(uint64_t(1), total_effective >> 21)); // / 2 MB
        size_t n = 1;
        while (n < raw) n <<= 1;           // next power of 2 >= raw
        n = std::max(n, size_t(16));
        n = std::min(n, std::min(max_parts, WRITER_CACHE_CAP));

        // For multi-file inputs where average file is small (< 50 MB uncompressed),
        // cap n at next_pow2(n_files): redundant genomes don't add unique k-mer content
        // proportionally, so total_effective overestimates partition load and leads to
        // too many tiny partitions with excessive phase2 overhead.
        const size_t n_files = cfg.input_files.size();
        if (n_files > 1) {
            const uint64_t avg_effective = total_effective / n_files;
            constexpr uint64_t SMALL_FILE_THRESHOLD = 50ULL << 20; // 50 MB
            if (avg_effective < SMALL_FILE_THRESHOLD) {
                size_t file_cap = 1;
                while (file_cap < n_files) file_cap <<= 1; // next_pow2(n_files)
                n = std::min(n, file_cap);
            }
        }
        cfg.num_partitions = static_cast<uint32_t>(n);
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
                  << "  l=" << cfg.l
                  << "  partitions=" << cfg.num_partitions
                  << "  threads=" << cfg.num_threads;
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
                fs::remove(cfg.work_dir + "hash_" + std::to_string(p) + ".superkmers");
        }
    }

    return rc;
}
