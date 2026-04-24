#include "CLI_prot.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

void print_prot_usage(const char* prog) {
    std::cerr <<
        "tuna-prot — protein k-mer counter\n\n"
        "Usage: " << prog << " [options] <file1> [file2 ...]\n\n"
        "Input:\n"
        "  Files may be plain FASTA or gzip FASTA (.gz).\n"
        "  Prefix a file with '@' to treat it as a file-of-filenames (one path per line).\n\n"
        "Options:\n"
        "  -k <int>   k-mer length              [default: 7,  range: 4-12]\n"
        "  -m <int>   minimizer length           [default: 5,  range: 2..k-1]\n"
        "  -n <int>   number of partitions       [default: auto, must be power of 2]\n"
        "  -t <int>   number of threads          [default: 4]\n"
        "  -o <file>  output file ('-' = stdout) [default: stdout]\n"
        "  -w <dir>   work directory             [default: .]\n"
        "  -ram <GB>  RAM budget (GB)            [default: auto-detect]\n"
        "  -hp        hide progress              [default: off]\n"
        "  -kt        keep temp partition files  [default: off]\n"
        "  -tp        stop after phase 1         [default: off]\n"
        "  -h         show this help\n\n"
        "Output: TSV with columns <kmer>\\t<count>, one k-mer per line.\n";
}

bool parse_prot_cli(int argc, char** argv, ProtConfig& cfg) {
    if (argc < 2) { print_prot_usage(argv[0]); return false; }

    for (int i = 1; i < argc; ++i) {
        const std::string a(argv[i]);

        if (a == "-h" || a == "--help") {
            print_prot_usage(argv[0]);
            return false;
        }

        auto need_next = [&](const char* flag) -> const char* {
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string(flag) + " requires an argument");
            return argv[++i];
        };

        if (a == "-k")   { cfg.k = static_cast<uint16_t>(std::atoi(need_next("-k"))); }
        else if (a == "-m")   { cfg.m = static_cast<uint16_t>(std::atoi(need_next("-m"))); }
        else if (a == "-n")   { cfg.num_partitions = static_cast<uint32_t>(std::atoi(need_next("-n"))); }
        else if (a == "-t")   { cfg.num_threads = static_cast<uint32_t>(std::atoi(need_next("-t"))); }
        else if (a == "-o")   { cfg.output_path = need_next("-o"); }
        else if (a == "-w")   { cfg.work_dir = need_next("-w"); }
        else if (a == "-ram") {
            const double gb = std::atof(need_next("-ram"));
            cfg.ram_budget_bytes = static_cast<uint64_t>(gb * (1ull << 30));
        }
        else if (a == "-hp")  { cfg.hide_progress = true; }
        else if (a == "-kt")  { cfg.keep_tmp = true; }
        else if (a == "-tp")  { cfg.stop_after_phase1 = true; }
        else if (a[0] != '-') {
            // Input file or @fof.
            if (!a.empty() && a[0] == '@') {
                // File-of-filenames.
                const std::string fof_path = a.substr(1);
                std::ifstream fof(fof_path);
                if (!fof) throw std::invalid_argument("Cannot open fof: " + fof_path);
                std::string line;
                while (std::getline(fof, line))
                    if (!line.empty()) cfg.input_files.push_back(line);
            } else {
                cfg.input_files.push_back(a);
            }
        } else {
            throw std::invalid_argument("Unknown option: " + a);
        }
    }

    if (cfg.input_files.empty())
        throw std::invalid_argument("No input files specified");
    if (cfg.k < 4 || cfg.k > 12)
        throw std::invalid_argument("-k must be in [4, 12]");
    if (cfg.m < 2 || cfg.m >= cfg.k)
        throw std::invalid_argument("-m must be in [2, k-1]");
    if (cfg.num_threads == 0) cfg.num_threads = 1;

    // Default work dir.
    if (cfg.work_dir.empty()) cfg.work_dir = ".";
    if (cfg.work_dir.back() != '/') cfg.work_dir += '/';

    return true;
}
