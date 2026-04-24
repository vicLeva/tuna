#include "CLI_prot.hpp"
#include "pipeline_prot.hpp"
#include "config_prot.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    ProtConfig cfg;
    try {
        if (!parse_prot_cli(argc, argv, cfg)) return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        print_prot_usage(argv[0]);
        return 1;
    }

    // Ensure work directory exists.
    try {
        std::filesystem::create_directories(cfg.work_dir);
    } catch (const std::exception& e) {
        std::cerr << "Error creating work dir " << cfg.work_dir << ": " << e.what() << "\n";
        return 1;
    }

    try {
        return prot::dispatch_run_prot(cfg);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
