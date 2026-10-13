#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

struct PartitionStats {
    uint64_t seqs = 0, kmers = 0, superkmers = 0;
    // Phase 1 started in memory, crossed the memory budget, and drained every
    // partition to disk. Phase 2 must then read from disk.
    bool spilled      = false;
    // A spill was needed but could not be completed (partition files could not
    // be opened or written). The run has to be reported as failed.
    bool spill_failed = false;
};

// The minimizer length tuna uses when -m is not given.
//
// Measured by scripts/bench_autom.sh over four regimes (500 bacterial
// assemblies, one human assembly, a read set and a metagenome) at nine values
// of k. Summed wall time picks 13 below k=31 and 17 at or above it; that rule
// is within 3.1% of the best m at every k measured and 0.6% on average. The
// per-k optima are flatter than they look -- at k=63 everything from m=15 to
// m=25 is within 1% -- so a two-step rule beats fitting each k, which would be
// fitting noise.
//
// Clamped to k-2 so small k stay valid: k=11 gives 9, k=13 gives 11.
constexpr uint16_t auto_m(uint16_t k) noexcept
{
    uint16_t m = (k < 31) ? uint16_t(13) : uint16_t(17);
    const uint16_t hi = (k > 2) ? uint16_t(k - 2) : uint16_t(1);
    if (m > hi) m = hi;
    if ((m & 1u) == 0u) --m;            // the dispatch table carries odd m only
    return m;
}

enum class SuperkmerDedupMode : uint8_t {
    Auto,
    On,
    Off
};

struct Config {
    std::vector<std::string> input_files;
    std::string  output_file;
    std::string  work_dir;            // directory for temp partition files
    uint16_t     k              = 31;
    uint16_t     l              = 0;  // 0 = auto: resolved to auto_m(k) by the CLI
    uint32_t     num_partitions = 0;    // 0 = auto-tune from input size (set in main)
    uint32_t     num_threads    = 1;  // worker threads: phase 1 over files, phase 2 over partitions
    uint32_t     ci             = 1;  // minimum count to report
    uint64_t     cx             = std::numeric_limits<uint64_t>::max(); // maximum count
    bool         hide_progress  = false;
    bool         keep_tmp      = false; // skip cleanup of partition files (useful for benchmarking)
    bool         partition_only = false; // exit after phase 1 (for benchmarking partition speed)
    bool         count_only    = false; // skip output writing after counting (for benchmarking)
    bool         debug_stats   = false; // print per-partition table stats + write minimizer coverage CSV
    bool         output_kff    = false; // write output in KFF binary format instead of TSV
    bool         canonical     = true;  // count canonical k-mers (lex-min of k-mer and RC); set false with -b
    uint64_t     ram_budget_bytes = 0; // 0 = auto-detect available RAM
    SuperkmerDedupMode dedup_mode = SuperkmerDedupMode::Auto;
};
