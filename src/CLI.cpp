#include "CLI.hpp"

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>


void print_usage(const char* prog)
{
    std::cerr <<
        "tuna — streaming k-mer counter\n"
        "\n"
        "Usage:\n"
        "  " << prog << " [options] <input1.fa [input2.fa ...]> <output_file>\n"
        "  " << prog << " [options] @<input_list_file>          <output_file>\n"
        "\n"
        "Options:\n"
        "  -k  <int>   k-mer length               [default: 31]\n"
        "              odd values in [11,31]  (k ≤ 32, fits in 64-bit word)\n"
        "  -m  <int>   minimizer length            [default: 21]\n"
        "              odd values in [9,k-2]  (must be odd and < k)\n"
        "              m=21 is a good default; use m=23-25 for highly repetitive\n"
        "              or low-complexity data (e.g. individual human genomes)\n"
        "  -n  <int>   number of partitions        [default: auto, ~2 MB input/partition]\n"
        "  -t  <int>   worker threads              [default: 1]\n"
        "              phase 1: parallel over input files (capped at #files)\n"
        "              phase 2: parallel over partitions  (capped at -n)\n"
        "  -ci <int>   minimum k-mer count         [default: 1]\n"
        "  -cx <int>   maximum k-mer count         [default: max]\n"
        "  -w  <dir>   working directory for temp files\n"
        "              [default: tuna_tmp/ next to output file]\n"
        "\n"
        "  -hp         hide progress messages\n"
        "  -kt         keep temp partition files after run (useful for benchmarking)\n"
        "  -tp         stop after partitioning (phase 1 only, for benchmarking)\n"
        "  -dbg        debug stats: per-partition table summary + minimizer coverage\n"
        "              CSV written to <work_dir>/debug_min_coverage.csv\n"
        "  -h/--help   show this help\n"
        "  -v/--version print version\n"
        "\n"
        "Output: plain text, one k-mer per line: <kmer>\\t<count>\n";
}


bool parse_args(int argc, char* argv[], Config& cfg)
{
    if (argc < 2) return false;

    std::vector<std::string> positionals;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") return false;

        if (arg == "-v" || arg == "--version") {
            std::cerr << "tuna " << TUNA_VERSION << "\n";
            std::exit(0);
        }

        // Consume the next token as the value for 'flag'.
        auto next_val = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "tuna: error: " << flag << " requires an argument\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "-k") {
            const char* v = next_val("-k"); if (!v) return false;
            cfg.k = static_cast<uint16_t>(std::stoul(v));
        } else if (arg == "-m") {
            const char* v = next_val("-m"); if (!v) return false;
            cfg.l = static_cast<uint16_t>(std::stoul(v));
        } else if (arg == "-n") {
            const char* v = next_val("-n"); if (!v) return false;
            cfg.num_partitions = static_cast<uint32_t>(std::stoul(v));
        } else if (arg == "-ci") {
            const char* v = next_val("-ci"); if (!v) return false;
            cfg.ci = static_cast<uint32_t>(std::stoul(v));
        } else if (arg == "-cx") {
            const char* v = next_val("-cx"); if (!v) return false;
            cfg.cx = std::stoull(v);
        } else if (arg == "-t") {
            const char* v = next_val("-t"); if (!v) return false;
            cfg.num_threads = static_cast<uint32_t>(std::stoul(v));
        } else if (arg == "-w") {
            const char* v = next_val("-w"); if (!v) return false;
            cfg.work_dir = v;
        } else if (arg == "-dbg") {
            cfg.debug_stats = true;
        } else if (arg == "-hp") {
            cfg.hide_progress = true;
        } else if (arg == "-kt") {
            cfg.keep_tmp = true;
        } else if (arg == "-tp") {
            cfg.partition_only = true;
        } else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "tuna: error: unknown option: " << arg << "\n";
            return false;
        } else {
            positionals.push_back(arg);
        }
    }

    if (positionals.size() < 2) {
        std::cerr << "tuna: error: need at least one input file and an output file\n";
        return false;
    }

    // Last positional is the output file; everything before is input.
    cfg.output_file = positionals.back();
    positionals.pop_back();

    // Single positional starting with '@' is a file listing input paths.
    if (positionals.size() == 1 && positionals[0][0] == '@') {
        std::string list_path = positionals[0].substr(1);
        // Expand a leading ~/ since the shell won't do it when ~ isn't the
        // first character of the token (e.g. @~/foo stays literal).
        if (list_path.size() >= 2 && list_path[0] == '~' && list_path[1] == '/') {
            const char* home = std::getenv("HOME");
            if (home) list_path = std::string(home) + list_path.substr(1);
        }
        std::ifstream list_f(list_path);
        if (!list_f) {
            std::cerr << "tuna: error: cannot open input list: " << list_path << "\n";
            return false;
        }
        std::string line;
        while (std::getline(list_f, line))
            if (!line.empty()) cfg.input_files.push_back(line);
    } else {
        cfg.input_files = std::move(positionals);
    }

    if (cfg.input_files.empty()) {
        std::cerr << "tuna: error: no input files\n";
        return false;
    }

    if (cfg.k % 2 == 0 || cfg.k < 11 || cfg.k > 31) {
        std::cerr << "tuna: error: k-mer length (-k " << cfg.k
                  << ") must be an odd value in [11, 31]\n";
        return false;
    }

    if (cfg.l % 2 == 0 || cfg.l < 9) {
        std::cerr << "tuna: error: minimizer length (-m " << cfg.l
                  << ") must be an odd value >= 9\n";
        return false;
    }

    if (cfg.l >= cfg.k) {
        std::cerr << "tuna: error: minimizer length (-m " << cfg.l
                  << ") must be strictly less than k-mer length (-k " << cfg.k << ")\n";
        return false;
    }

    if (cfg.num_threads == 0) {
        std::cerr << "tuna: error: number of threads (-t) must be > 0\n";
        return false;
    }

    return true;
}
