#pragma once

// Phase 2+3 — partition files → counted k-mers → output.
//
// Backend: ankerl::unordered_dense::map (flat Robin Hood hash map).
// K-mers extracted as canonical uint64_t (k≤32) or __uint128_t (k≤64),
// with fwd/rc sliding window: O(k) init per superkmer, O(1) per base.

#include "Config.hpp"
#include "kff_output.hpp"
#include "superkmer_io.hpp"

#include <array>
#include <fstream>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <type_traits>
#include <vector>
#include <exception>
#include <stdexcept>
#include <charconv>
#include <iostream>

#include "ankerl/unordered_dense.h"


// Returns the smallest power of two >= v (returns 1 for v == 0).
inline size_t next_pow2(size_t v) noexcept {
    size_t p = 1;
    while (p < v) p <<= 1;
    return p;
}


// ─── K-mer types ──────────────────────────────────────────────────────────────

template <uint16_t k>
using kmer_t = std::conditional_t<(k <= 32), uint64_t, __uint128_t>;

struct kmer128_hash {
    using is_avalanching = void;
    auto operator()(__uint128_t x) const noexcept -> uint64_t {
        const uint64_t lo = static_cast<uint64_t>(x);
        const uint64_t hi = static_cast<uint64_t>(x >> 64);
        return ankerl::unordered_dense::hash<uint64_t>{}(lo ^ (hi * UINT64_C(0x9e3779b97f4a7c15)));
    }
};

template <uint16_t k>
using table_hash_t = std::conditional_t<(k <= 32),
    ankerl::unordered_dense::hash<uint64_t>,
    kmer128_hash>;

template <uint16_t k>
using table_t = ankerl::unordered_dense::map<kmer_t<k>, uint32_t, table_hash_t<k>>;


// ─── Debug stats ──────────────────────────────────────────────────────────────

struct PartitionDebugInfo {
    uint64_t init_sz    = 0;
    uint64_t n_inserted = 0;
    uint64_t n_unique   = 0;
    uint64_t table_cap  = 0;   // bucket_count (index layer size)
    uint64_t n_buckets  = 0;
};


// ─── K-mer serialization helpers ─────────────────────────────────────────────

static const char KMER_BASES[] = "ACGT";

template <uint16_t k>
void kmer_to_str(kmer_t<k> kmer, std::string& dst)
{
    dst.resize(k);
    for (int i = static_cast<int>(k) - 1; i >= 0; --i) {
        dst[static_cast<size_t>(i)] = KMER_BASES[static_cast<uint8_t>(kmer & 3u)];
        kmer >>= 2;
    }
}

// Write k-mer as 2-bit MSB-first packed bytes (KFF format).
// KMER_BYTES = ceil(k/4); padding bits (if any) at LSB of last byte.
template <uint16_t k>
void kmer_to_packed(kmer_t<k> kmer, uint8_t* dst)
{
    constexpr size_t KMER_BYTES = (k + 3u) / 4u;
    constexpr int    pad        = static_cast<int>(KMER_BYTES * 8u) - 2 * static_cast<int>(k);
    auto aligned = kmer << pad;
    for (size_t i = 0; i < KMER_BYTES; ++i)
        dst[i] = static_cast<uint8_t>(aligned >> (8u * (KMER_BYTES - 1u - i)));
}


// ─── Counting brick ───────────────────────────────────────────────────────────
//
// Extracts canonical k-mers from packed superkmers into an unordered_dense map.
// Sliding window: fwd = (fwd << 2) | b;  rc = (rc >> 2) | (complement(b) << 2*(k-1))
// canonical = min(fwd, rc).  O(k) init per superkmer, O(1) per additional base.

template <uint16_t k, uint16_t m, typename Reader = SuperkmerReader<k, m>>
uint64_t count_partition(Reader& reader, table_t<k>& table)
{
    using kmer_type = kmer_t<k>;

    static constexpr kmer_type kmer_mask = (2u * k < sizeof(kmer_type) * 8u)
        ? (kmer_type(1) << (2u * k)) - kmer_type(1)
        : ~kmer_type(0);
    static constexpr unsigned rc_shift = 2u * (k - 1u);

    uint64_t inserted = 0;

    while (reader.next()) {
        const uint8_t* packed = reader.packed_data();
        const size_t   len    = reader.size();
        if (len < k) continue;

        // Init fwd and rc from the first k bases.
        kmer_type fwd = 0, rc = 0;
        for (size_t i = 0; i < k; ++i) {
            const uint8_t b = (packed[i >> 2] >> (6u - 2u * (i & 3u))) & 3u;
            fwd = (fwd << 2) | b;
            rc  = (rc >> 2) | (kmer_type(b ^ 3u) << rc_shift);
        }
        ++table[std::min(fwd, rc)];
        ++inserted;

        // Slide through the remaining bases.
        const uint8_t* byte_ptr = packed + (k >> 2);
        int shift = static_cast<int>(6u - 2u * (k & 3u));
        for (size_t i = k; i < len; ++i) {
            const uint8_t b = (*byte_ptr >> shift) & 3u;
            shift -= 2;
            if (shift < 0) { shift = 6; ++byte_ptr; }
            fwd = ((fwd << 2) & kmer_mask) | b;
            rc  = (rc >> 2) | (kmer_type(b ^ 3u) << rc_shift);
            ++table[std::min(fwd, rc)];
            ++inserted;
        }
    }
    return inserted;
}


// ─── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
uint64_t write_counts(
    table_t<k>&     table,
    const Config&   cfg,
    std::string&    chunk,
    std::ofstream&  out,
    std::mutex&     out_mutex)
{
    constexpr size_t WRITE_BATCH = 1u << 20;
    if (chunk.capacity() < WRITE_BATCH + 64) chunk.reserve(WRITE_BATCH + 64);

    const auto flush_chunk = [&]() {
        std::lock_guard<std::mutex> g(out_mutex);
        out.write(chunk.data(), chunk.size());
        chunk.clear();
    };

    std::string label;
    uint64_t written = 0;

    for (const auto& [key, cnt] : table) {
        if (cnt < cfg.ci || cnt > cfg.cx) continue;
        kmer_to_str<k>(key, label);
        chunk += label;
        chunk += '\t';
        char cnt_buf[32];
        const auto [ptr, ec] = std::to_chars(std::begin(cnt_buf), std::end(cnt_buf), cnt);
        if (ec == std::errc()) chunk.append(cnt_buf, ptr);
        else chunk += std::to_string(cnt);
        chunk += '\n';
        ++written;
        if (chunk.size() >= WRITE_BATCH) flush_chunk();
    }
    if (!chunk.empty()) flush_chunk();
    return written;
}


// ─── KFF output brick ─────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
uint64_t write_counts_kff(
    table_t<k>&   table,
    const Config& cfg,
    KffOutput&    kff_out)
{
    constexpr size_t KMER_BYTES  = (k + 3u) / 4u;
    constexpr size_t BATCH_KMERS = (size_t(1) << 20) / (KMER_BYTES + 4u);

    std::vector<uint8_t>  seq_buf;
    std::vector<uint32_t> cnt_buf;
    seq_buf.reserve(BATCH_KMERS * KMER_BYTES);
    cnt_buf.reserve(BATCH_KMERS);

    const auto flush = [&]() {
        if (!cnt_buf.empty()) {
            kff_out.write_batch(seq_buf.data(), cnt_buf.data(), cnt_buf.size());
            seq_buf.clear();
            cnt_buf.clear();
        }
    };

    uint64_t written = 0;
    for (const auto& [key, cnt] : table) {
        if (cnt < cfg.ci || cnt > cfg.cx) continue;
        seq_buf.resize(seq_buf.size() + KMER_BYTES);
        kmer_to_packed<k>(key, seq_buf.data() + seq_buf.size() - KMER_BYTES);
        cnt_buf.push_back(cnt);
        ++written;
        if (cnt_buf.size() >= BATCH_KMERS) flush();
    }
    flush();
    return written;
}


// ─── Callback output brick ────────────────────────────────────────────────────

template <uint16_t k, uint16_t m, typename Callback>
uint64_t write_counts_callback(
    table_t<k>&   table,
    const Config& cfg,
    Callback&     cb)
{
    uint64_t written = 0;
    std::string label;

    for (const auto& [key, cnt] : table) {
        if (cnt < cfg.ci || cnt > cfg.cx) continue;
        if constexpr (std::is_invocable_v<Callback, kmer_t<k>, uint32_t>) {
            cb(key, cnt);
        } else {
            kmer_to_str<k>(key, label);
            cb(std::string_view(label), cnt);
        }
        ++written;
    }
    return written;
}


// ─── Debug output helper ──────────────────────────────────────────────────────

inline void emit_debug_stats(
    const std::vector<PartitionDebugInfo>& part_infos,
    size_t n_parts,
    const Config& cfg)
{
    std::cerr << "\n[debug] per-partition table stats (first 20):\n";
    std::cerr << "  part   init_sz   table_cap   n_inserted    n_unique  load%\n";
    for (size_t p = 0; p < std::min(n_parts, size_t(20)); ++p) {
        const auto& d = part_infos[p];
        const double lf = d.table_cap ? 100.0 * d.n_unique / d.table_cap : 0.0;
        std::cerr << "  " << std::setw(5)  << p
                  << "  " << std::setw(8)  << d.init_sz
                  << "  " << std::setw(9)  << d.table_cap
                  << "  " << std::setw(11) << d.n_inserted
                  << "  " << std::setw(10) << d.n_unique
                  << "  " << std::fixed << std::setprecision(1) << std::setw(5) << lf << "%\n";
    }
    if (n_parts > 20) std::cerr << "  ... (" << (n_parts - 20) << " more partitions)\n";

    double   sum_lf = 0.0, min_lf = 2.0, max_lf = -1.0;
    uint64_t sum_unique = 0, min_unique = uint64_t(-1), max_unique = 0;
    uint64_t sum_init = 0,   min_init   = uint64_t(-1), max_init   = 0;
    size_t   n_valid = 0;

    for (size_t p = 0; p < n_parts; ++p) {
        const auto& d = part_infos[p];
        if (d.table_cap == 0) continue;
        ++n_valid;
        const double lf = static_cast<double>(d.n_unique) / d.table_cap;
        sum_lf    += lf;
        min_lf     = std::min(min_lf, lf);
        max_lf     = std::max(max_lf, lf);
        sum_init  += d.init_sz;
        min_init   = std::min(min_init, d.init_sz);
        max_init   = std::max(max_init, d.init_sz);
        sum_unique += d.n_unique;
        min_unique  = std::min(min_unique, d.n_unique);
        max_unique  = std::max(max_unique, d.n_unique);
    }

    if (n_valid > 0) {
        const double mean_lf     = sum_lf / n_valid;
        const double mean_unique = static_cast<double>(sum_unique) / n_valid;
        const double mean_init   = static_cast<double>(sum_init) / n_valid;
        std::cerr << "\n[debug] aggregate table stats (" << n_parts << " partitions):\n"
                  << std::fixed;
        std::cerr << "  load_factor:   mean=" << std::setprecision(3) << mean_lf
                  << "  min=" << min_lf << "  max=" << max_lf << "\n";
        std::cerr << "  n_unique/part: mean=" << std::setprecision(0) << mean_unique
                  << "  min=" << min_unique << "  max=" << max_unique << "\n";
        std::cerr << "  init_sz/part:  mean=" << mean_init
                  << "  min=" << min_init << "  max=" << max_init << "\n";
        std::cerr << "dbg_load_mean: "   << std::setprecision(6) << mean_lf    << "\n"
                  << "dbg_load_min: "    << min_lf                              << "\n"
                  << "dbg_load_max: "    << max_lf                              << "\n"
                  << "dbg_unique_mean: " << std::setprecision(0) << mean_unique << "\n";
    }

    const std::string csv_path = cfg.work_dir + "debug_table_stats.csv";
    std::ofstream csv(csv_path);
    if (csv) {
        csv << "partition_id,init_sz,table_cap,n_inserted,n_unique,load_factor\n";
        for (size_t p = 0; p < n_parts; ++p) {
            const auto& d = part_infos[p];
            const double lf = d.table_cap
                ? static_cast<double>(d.n_unique) / d.table_cap : 0.0;
            csv << p << "," << d.init_sz << "," << d.table_cap << ","
                << d.n_inserted << "," << d.n_unique << ","
                << std::fixed << std::setprecision(6) << lf << "\n";
        }
        std::cerr << "[debug] table stats CSV: " << csv_path << "\n";
    }
}


// ─── init_sz helper ───────────────────────────────────────────────────────────

inline size_t compute_init_sz(uint64_t cal, uint64_t total_kmers, size_t n_parts) noexcept
{
    if (cal > 0)
        return std::clamp(static_cast<size_t>(cal),
                          size_t(1u << 15), size_t(1u << 22));
    const size_t per_part = total_kmers > 0
        ? static_cast<size_t>(total_kmers / n_parts)
        : size_t(1u << 27) / n_parts;
    return std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
}


// ─── Callback counting harnesses ──────────────────────────────────────────────

template <uint16_t k, uint16_t m, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    Callback&&                cb)
{
    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t>   next_part{0};
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t) {
        try {
            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const size_t init_sz = compute_init_sz(
                    calibrated_unique.load(std::memory_order_relaxed), total_kmers, n_parts);
                table_t<k> table;
                table.reserve(init_sz);

                uint64_t ins;
                {
                    MemoryReader<k, m> reader(part_bufs[p]);
                    ins = count_partition<k, m, MemoryReader<k, m>>(reader, table);
                }
                { std::string tmp; part_bufs[p].swap(tmp); }
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                const uint64_t unique = static_cast<uint64_t>(table.size());
                if (unique > 0) {
                    uint64_t expected = 0;
                    calibrated_unique.compare_exchange_strong(
                        expected, unique, std::memory_order_relaxed, std::memory_order_relaxed);
                }

                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);
                total_written.fetch_add(wrt, std::memory_order_relaxed);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(worker_error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);

    return { total_inserted.load(), total_written.load() };
}


template <uint16_t k, uint16_t m, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback(
    const Config& cfg,
    uint64_t      total_kmers,
    Callback&&    cb)
{
    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t>   next_part{0};
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t) {
        try {
            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const std::string path = partition_path(cfg.work_dir, p);
                SuperkmerReader<k, m> reader(path);
                if (!reader.ok())
                    throw std::runtime_error("tuna: cannot open partition file: " + path);

                const size_t init_sz = compute_init_sz(
                    calibrated_unique.load(std::memory_order_relaxed), total_kmers, n_parts);
                table_t<k> table;
                table.reserve(init_sz);

                const uint64_t ins = count_partition<k, m>(reader, table);
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                const uint64_t unique = static_cast<uint64_t>(table.size());
                if (unique > 0) {
                    uint64_t expected = 0;
                    calibrated_unique.compare_exchange_strong(
                        expected, unique, std::memory_order_relaxed, std::memory_order_relaxed);
                }
                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);
                total_written.fetch_add(wrt, std::memory_order_relaxed);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(worker_error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);

    return { total_inserted.load(), total_written.load() };
}


// ─── count_and_write (disk mode) ─────────────────────────────────────────────

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    uint64_t       total_kmers,
    std::ofstream* out,
    KffOutput*     kff_out)
{
    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t>   next_part{0};
    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t) {
        try {
            std::string chunk;
            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const std::string path = partition_path(cfg.work_dir, p);
                SuperkmerReader<k, m> reader(path);
                if (!reader.ok())
                    throw std::runtime_error("tuna: cannot open partition file: " + path);

                const size_t init_sz = compute_init_sz(
                    calibrated_unique.load(std::memory_order_relaxed), total_kmers, n_parts);
                table_t<k> table;
                table.reserve(init_sz);

                const uint64_t ins = count_partition<k, m>(reader, table);
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                const uint64_t unique = static_cast<uint64_t>(table.size());
                if (unique > 0) {
                    uint64_t expected = 0;
                    calibrated_unique.compare_exchange_strong(
                        expected, unique, std::memory_order_relaxed, std::memory_order_relaxed);
                }

                if (cfg.debug_stats) {
                    auto& d = part_infos[p];
                    d.init_sz    = init_sz;
                    d.n_inserted = ins;
                    d.n_unique   = unique;
                    d.table_cap  = table.bucket_count();
                    d.n_buckets  = table.bucket_count();
                }

                const uint64_t wrt = kff_out
                    ? write_counts_kff<k, m>(table, cfg, *kff_out)
                    : write_counts<k, m>(table, cfg, chunk, *out, out_mutex);
                total_written.fetch_add(wrt, std::memory_order_relaxed);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(worker_error_mutex);
            if (!worker_error) worker_error = std::current_exception();
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);

    if (cfg.debug_stats) emit_debug_stats(part_infos, n_parts, cfg);

    return { total_inserted.load(), total_written.load() };
}


// ─── count_and_write_mem (in-memory mode) ────────────────────────────────────

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> count_and_write_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    std::ofstream*            out,
    KffOutput*                kff_out)
{
    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t>   next_part{0};
    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t) {
        std::string chunk;
        while (true) {
            const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
            if (p >= n_parts) break;

            const size_t init_sz = compute_init_sz(
                calibrated_unique.load(std::memory_order_relaxed), total_kmers, n_parts);
            table_t<k> table;
            table.reserve(init_sz);

            uint64_t ins;
            {
                MemoryReader<k, m> reader(part_bufs[p]);
                ins = count_partition<k, m, MemoryReader<k, m>>(reader, table);
            }
            { std::string tmp; part_bufs[p].swap(tmp); }
            total_inserted.fetch_add(ins, std::memory_order_relaxed);

            const uint64_t unique = static_cast<uint64_t>(table.size());
            if (unique > 0) {
                uint64_t expected = 0;
                calibrated_unique.compare_exchange_strong(
                    expected, unique, std::memory_order_relaxed, std::memory_order_relaxed);
            }

            if (cfg.debug_stats) {
                auto& d = part_infos[p];
                d.init_sz    = init_sz;
                d.n_inserted = ins;
                d.n_unique   = unique;
                d.table_cap  = table.bucket_count();
                d.n_buckets  = table.bucket_count();
            }

            const uint64_t wrt = kff_out
                ? write_counts_kff<k, m>(table, cfg, *kff_out)
                : write_counts<k, m>(table, cfg, chunk, *out, out_mutex);
            total_written.fetch_add(wrt, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    if (cfg.debug_stats) emit_debug_stats(part_infos, n_parts, cfg);

    return { total_inserted.load(), total_written.load() };
}
