#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

enum class PartitionStrategy {
    HASH,      // simple minimizer-hash % num_partitions — no pre-scan (default)
    KMTRICKS,  // kmtricks-style: rolling minimizer hash + load-balanced table
    KMC,       // KMC norm-filtered signature + load-balanced table
};

struct Config {
    std::vector<std::string> input_files;
    std::string  output_file;
    std::string  work_dir;            // directory for temp partition files
    uint16_t     k              = 31;
    uint16_t     l              = 17;
    uint32_t     num_partitions = 32;
    uint32_t     num_threads    = 1;  // worker threads: phase 1 over files, phase 2 over partitions
    uint32_t     ci             = 1;  // minimum count to report
    uint64_t     cx             = std::numeric_limits<uint64_t>::max(); // maximum count
    bool         hide_progress  = false;
    PartitionStrategy strategy  = PartitionStrategy::HASH;
    uint16_t     kmc_sig_len   = 9;    // KMC signature length (5–11, default 9)
    bool         keep_tmp      = false; // skip cleanup of partition files (useful for benchmarking)
    bool         partition_only = false; // exit after phase 1 (for benchmarking partition speed)

    // Returns the strategy prefix used in partition file names.
    std::string partition_prefix() const {
        switch (strategy) {
            case PartitionStrategy::KMC:      return "kmc";
            case PartitionStrategy::KMTRICKS: return "kmtricks";
            case PartitionStrategy::HASH:     return "hash";
        }
        return "hash";
    }
};
