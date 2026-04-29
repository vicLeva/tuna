#pragma once

// tuna — public C++ API
//
// Count canonical k-mers in FASTA/FASTQ files (plain or gzipped).
//
// ── Quick start ───────────────────────────────────────────────────────────────
//
//   // Collect into an unordered_map (safe, single-threaded container):
//   auto kmers = tuna::count_to<31>({"genome.fa"});
//   // kmers is std::unordered_map<std::string, uint32_t>
//
//   // Streaming callback (may be called from multiple threads):
//   tuna::count<31>({"genome.fa"}, [](std::string_view kmer, uint32_t cnt) {
//       process(kmer, cnt);
//   });
//
// ── CMake integration ─────────────────────────────────────────────────────────
//
//   add_subdirectory(tuna)
//   target_link_libraries(my_target PRIVATE tuna::tuna)
//
// ── Template parameters ───────────────────────────────────────────────────────
//
//   k — k-mer length in [2, 256]  (fixed at compile time via FIXED_K)
//       Kmer<k> uses ceil(k/32) uint64 words (1 for k≤32, 2 for k≤64, …).
//   m — minimizer length in [1, k-1]; default 21
//       Odd values recommended; use m=21–25 for standard genomics work.

#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <filesystem>

// Internal implementation headers (all template code, must be compiled inline).
#include "pipeline.hpp"
#include "setup.hpp"
#include "superkmer_io.hpp"   // SuperkmerWriter (for sizeof in auto-tune)

namespace tuna {

// ─── K-mer type ───────────────────────────────────────────────────────────────
//
// tuna::Kmer<k> is the packed binary representation of a canonical k-mer.
// Use it as the argument type in binary-mode callbacks (see count() below).
// uint64_t for k ≤ 32, __uint128_t for k ≤ 64.  2 bits per base, MSB-first.

template <uint16_t k>
using Kmer = kmer_t<k>;

// ─── Options ──────────────────────────────────────────────────────────────────

struct Options {
    uint32_t    threads    = 1;     // worker threads
    uint32_t    min_count  = 1;     // only yield k-mers with count ≥ min_count
    uint32_t    max_count  = std::numeric_limits<uint32_t>::max(); // ...and ≤ max_count
    uint32_t    partitions = 0;     // 0 = auto-tune (recommended)
    uint32_t    ram_gb     = 0;     // RAM budget in GB; 0 = auto-detect available RAM
    std::string work_dir;           // "" = auto-managed temp dir, deleted on return
};


// ─── count() ──────────────────────────────────────────────────────────────────
//
// Invokes cb(kmer, count) for each canonical k-mer that passes the count filters.
//
// Two callback signatures are supported — the right path is selected at compile
// time based on what the callback accepts:
//
//   ASCII  (default): cb(std::string_view kmer, uint32_t count)
//   Binary           : cb(tuna::Kmer<k> kmer, uint32_t count)
//
// Binary mode skips the packed→string conversion and is slightly faster.
// Use kmer_to_str<k>(kmer, str) to decode a binary k-mer to an ASCII string.
//
// cb is called concurrently: with num_threads > 1, multiple worker threads process
// different partitions in parallel and each calls cb independently.  K-mers are
// disjoint across partitions, so no k-mer is reported twice.  The callback must be
// thread-safe (use a mutex if writing to shared state — see count_to for a
// ready-made example).
//
// Throws std::invalid_argument if files is empty.
// Throws std::runtime_error if the temp directory cannot be created.

template <uint16_t k, uint16_t m = 21, typename Callback>
void count(const std::vector<std::string>& files, Callback&& cb, Options opts = {})
{
    static_assert(k >= 2 && k <= 256, "tuna: k must be in [2, 256]");
    static_assert(m >= 1 && m < k,   "tuna: m must be in [1, k-1]");

    if (files.empty())
        throw std::invalid_argument("tuna::count: no input files");

    // ── Build internal config ──────────────────────────────────────────────
    Config cfg;
    cfg.input_files    = files;
    cfg.output_file    = "";
    cfg.k              = k;
    cfg.l              = m;
    cfg.num_threads       = opts.threads > 0 ? opts.threads : 1;
    cfg.ci                = opts.min_count;
    cfg.cx                = opts.max_count;
    cfg.num_partitions    = opts.partitions;
    cfg.ram_budget_bytes  = opts.ram_gb > 0
        ? static_cast<uint64_t>(opts.ram_gb) << 30 : 0;
    cfg.hide_progress     = true;

    static_assert(sizeof(SuperkmerWriter<k, m>) >= 8 && sizeof(SuperkmerWriter<k, m>) <= 64,
                  "SuperkmerWriter size out of expected range");
    if (cfg.num_partitions == 0)
        cfg.num_partitions = auto_tune_partitions(
            cfg.input_files, sizeof(SuperkmerWriter<k, m>), cfg.k, cfg.num_threads);
    cfg.num_partitions = round_pow2(cfg.num_partitions);

    // ── Work directory ─────────────────────────────────────────────────────
    namespace fs = std::filesystem;
    const bool own_work_dir = opts.work_dir.empty();
    if (own_work_dir) {
        // mkdtemp requires a writable char buffer ending in XXXXXX.
        std::string tmpl = (fs::temp_directory_path() / "tuna_XXXXXX").string();
        if (!mkdtemp(tmpl.data()))
            throw std::runtime_error(
                std::string("tuna: failed to create temp directory: ") + strerror(errno));
        cfg.work_dir = std::move(tmpl);
    } else {
        cfg.work_dir = opts.work_dir;
        fs::create_directories(cfg.work_dir);
    }
    if (cfg.work_dir.back() != '/') cfg.work_dir += '/';

    // ── Run pipeline and clean up ──────────────────────────────────────────
    try {
        run_callback<k, m>(cfg, std::forward<Callback>(cb));
    } catch (...) {
        if (own_work_dir) fs::remove_all(cfg.work_dir);
        throw;
    }
    if (own_work_dir) fs::remove_all(cfg.work_dir);
}


// ─── count_to() ───────────────────────────────────────────────────────────────
//
// Collects k-mer → count pairs into Container and returns it.
//
// Container must support emplace(std::string, uint32_t) or equivalent insertion.
// Defaults to std::unordered_map<std::string, uint32_t>.
//
// Thread safety: a mutex is used internally so Container does not need to be
// thread-safe.  For custom containers, emplace must accept (std::string, uint32_t).

template <uint16_t k, uint16_t m = 21,
          typename Container = std::unordered_map<std::string, uint32_t>>
Container count_to(const std::vector<std::string>& files, Options opts = {})
{
    Container result;
    std::mutex mtx;
    count<k, m>(
        files,
        [&](std::string_view kmer, uint32_t cnt) {
            std::lock_guard<std::mutex> lg(mtx);
            result.emplace(std::string(kmer), cnt);
        },
        opts);
    return result;
}

// ─── count_to_raw() ───────────────────────────────────────────────────────────
//
// Binary variant of count_to(): collects (Kmer<k>, count) pairs into a
// std::vector and returns it.  No packed→string conversion is performed.
//
// Use kmer_to_str<k>(kmer, str) to decode a k-mer back to ASCII if needed.

template <uint16_t k, uint16_t m = 21>
std::vector<std::pair<Kmer<k>, uint32_t>>
count_to_raw(const std::vector<std::string>& files, Options opts = {})
{
    std::vector<std::pair<Kmer<k>, uint32_t>> result;
    std::mutex mtx;
    count<k, m>(
        files,
        [&](Kmer<k> kmer, uint32_t cnt) {
            std::lock_guard<std::mutex> lg(mtx);
            result.emplace_back(kmer, cnt);
        },
        opts);
    return result;
}

} // namespace tuna
