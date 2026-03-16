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

#include <filesystem>
#include <iostream>
#include <limits>
#include <sys/resource.h>


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
    // Target: ~2 MB raw input per partition.  For compressed inputs file_size()
    // returns compressed bytes — an underestimate — so n is slightly inflated,
    // which is safe (smaller tables, better L3 fit).
    // n is rounded to the next power of 2 and clamped to [16, fd_budget] where
    // fd_budget = current RLIMIT_NOFILE - 32 (headroom for other FDs).
    if (cfg.num_partitions == 0) {
        struct rlimit rl{};
        size_t fd_limit = 1024;
        if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
            fd_limit = static_cast<size_t>(rl.rlim_cur);
        const size_t max_parts = (fd_limit > 32) ? fd_limit - 32 : 16;

        uint64_t total_bytes = 0;
        for (const auto& f : cfg.input_files) {
            std::error_code ec;
            total_bytes += fs::file_size(f, ec);
        }
        const size_t raw = static_cast<size_t>(
            std::max(uint64_t(1), total_bytes >> 21)); // / 2 MB
        size_t n = 1;
        while (n < raw) n <<= 1;           // next power of 2 >= raw
        n = std::max(n, size_t(16));
        n = std::min(n, max_parts);
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
