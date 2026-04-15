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
#include "setup.hpp"
#include "superkmer_io.hpp"

#include <filesystem>
#include <iostream>
#include <limits>
#include <fstream>
#include <sys/resource.h>
#include <exception>


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

    // SuperkmerWriter struct layout is independent of k/m (hdr_t only affects
    // serialized bytes, not struct fields) — any instantiation gives the same sizeof.
    static_assert(sizeof(SuperkmerWriter<31, 21>) >= 8 && sizeof(SuperkmerWriter<31, 21>) <= 64,
                  "SuperkmerWriter size out of expected range");

    for (const auto& input : cfg.input_files) {
        std::ifstream f(input, std::ios::binary);
        if (!f) {
            std::cerr << "tuna: error: cannot open input file: " << input << "\n";
            return 1;
        }
    }

    if (cfg.num_partitions == 0)
        cfg.num_partitions = auto_tune_partitions(
            cfg.input_files, sizeof(SuperkmerWriter<31, 21>), cfg.k, cfg.num_threads);

    // Round up to next power of 2 so partition_fn can use bitmask instead of division.
    cfg.num_partitions = round_pow2(cfg.num_partitions);

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

    int rc = 1;
    try {
        rc = dispatch(cfg.k, cfg.l, cfg);
    } catch (const std::exception& e) {
        std::cerr << "tuna: error: " << e.what() << "\n";
        rc = 1;
    } catch (...) {
        std::cerr << "tuna: error: unexpected failure\n";
        rc = 1;
    }

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
