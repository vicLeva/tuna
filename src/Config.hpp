#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <limits>

struct PartitionStats { uint64_t seqs = 0, kmers = 0, superkmers = 0; };

struct Config {
    std::vector<std::string> input_files;
    std::string  output_file;
    std::string  work_dir;            // directory for temp partition files
    uint16_t     k              = 31;
    uint16_t     l              = 21;
    uint32_t     num_partitions = 0;    // 0 = auto-tune from input size (set in main)
    uint32_t     num_threads    = 1;  // worker threads: phase 1 over files, phase 2 over partitions
    uint32_t     ci             = 1;  // minimum count to report
    uint64_t     cx             = std::numeric_limits<uint64_t>::max(); // maximum count
    bool         hide_progress  = false;
    bool         keep_tmp      = false; // skip cleanup of partition files (useful for benchmarking)
    bool         partition_only = false; // exit after phase 1 (for benchmarking partition speed)
    bool         phase2_only   = false; // skip phase 1 and count existing partition files in work_dir
    bool         count_only    = false; // count k-mers but skip output writing
    bool         debug_stats   = false; // print per-partition table stats + write minimizer coverage CSV
    bool         output_kff    = false; // write output in KFF binary format instead of TSV
    bool         lz4_buckets   = false; // compress disk-mode phase-1 bucket files with LZ4
    bool         lz4_shards    = false; // compress recursive phase-2 dedup shard files with LZ4
    bool         phase1_adaptive = false; // use experimental adaptive phase-1 scheduler for gz FASTQ
    uint64_t     ram_budget_bytes = 0; // 0 = auto-detect available RAM
};
