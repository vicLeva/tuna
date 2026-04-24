#pragma once
#include <cstdint>
#include <string>
#include <vector>

struct ProtConfig {
    std::vector<std::string> input_files;
    std::string              output_path;    // "-" or "" → stdout
    std::string              work_dir;

    uint16_t k            = 7;
    uint16_t m            = 5;
    uint32_t num_partitions = 0;  // 0 = auto
    uint32_t num_threads  = 4;

    uint64_t ram_budget_bytes = 0;  // 0 = detect

    bool hide_progress    = false;
    bool keep_tmp         = false;
    bool stop_after_phase1 = false;  // -tp flag
};
