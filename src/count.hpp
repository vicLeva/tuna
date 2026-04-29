#pragma once

// Phase 2+3 — partition files → counted k-mers → output.
//
// Three independently modifiable bricks:
//
//   count_partition<k,m>   — read one partition file, fill the iceberg hash table.
//   write_counts<k,m>      — drain one table to the output stream.
//   count_and_write<k,m>   — parallel harness (threading / scheduling).

#include "Config.hpp"
#include "kff_output.hpp"
#include "superkmer_io.hpp"

#include <array>
#include <fstream>
#include <iostream>
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
#include <bit>

#include "iceberg_table.h"


// Returns the smallest power of two >= v (returns 1 for v == 0).
inline size_t next_pow2(size_t v) noexcept {
    size_t p = 1;
    while (p < v) p <<= 1;
    return p;
}

// Compute log_slots for iceberg_init from an element-count estimate.
// Targets ~50% lv1 occupancy: nslots = next_pow2(2 * n).
// Minimum log_slots = 7 (SLOT_BITS+1 = 128 slots, 2 blocks).
inline uint64_t iceberg_log_slots(size_t n) noexcept {
    const size_t nslots = std::max(size_t(128), next_pow2(2 * n));
    return static_cast<uint64_t>(std::countr_zero(nslots));
}


// ─── K-mer type machinery ─────────────────────────────────────────────────────
//
// iceberg KeyType is uint64_t, so only k <= 32 is supported.

template <uint16_t k>
using kmer_t = uint64_t;

static_assert(true); // placeholder; actual k check is in count_partition

template <uint16_t k>
inline uint64_t kmer_mask_v() noexcept {
    if constexpr (k < 32)
        return (uint64_t(1) << (2 * k)) - 1;
    else
        return ~uint64_t(0);
}


// ─── Debug stats ──────────────────────────────────────────────────────────────

struct PartitionDebugInfo {
    uint64_t init_sz    = 0;
    uint64_t n_inserted = 0;
    uint64_t n_unique   = 0;
    uint64_t table_cap  = 0;   // nslots = 1 << log_slots
};


// ─── K-mer string conversion ──────────────────────────────────────────────────

template <uint16_t k>
inline void kmer_to_str(uint64_t key, std::string& out) noexcept {
    static constexpr char b[] = "ACGT";
    out.resize(k);
    for (int i = static_cast<int>(k) - 1; i >= 0; --i, key >>= 2)
        out[static_cast<size_t>(i)] = b[key & 3];
}

template <uint16_t k>
inline void kmer_to_packed_msb(uint64_t key, uint8_t* dst) noexcept {
    constexpr size_t KMER_BYTES = (k + 3) / 4;
    const unsigned shift = static_cast<unsigned>(KMER_BYTES * 8 - 2 * k);
    const uint64_t lj = key << shift;
    for (size_t i = 0; i < KMER_BYTES; ++i)
        dst[i] = static_cast<uint8_t>(lj >> (8 * (KMER_BYTES - 1 - i)));
}


// ─── Counting brick ───────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m, typename Reader = SuperkmerReader<k, m>>
uint64_t count_partition(
    Reader&         reader,
    iceberg_table&  table,
    PartitionDebugInfo* dbg = nullptr)
{
    static_assert(k <= 32, "iceberg backend: k must be <= 32 (KeyType is uint64_t)");

    const uint64_t mask = kmer_mask_v<k>();
    uint64_t inserted = 0;

    while (reader.next()) {
        const uint8_t* packed = reader.packed_data();
        const size_t   len    = reader.size();
        if (len < k) continue;

        // ── Init first k-mer ─────────────────────────────────────────────────
        uint64_t fwd = 0;
        for (uint16_t i = 0; i < k; ++i)
            fwd = (fwd << 2) | uint64_t((packed[i >> 2] >> (6 - 2 * (i & 3))) & 3u);

        uint64_t rc = 0;
        {
            uint64_t tmp = fwd;
            for (uint16_t i = 0; i < k; ++i, tmp >>= 2)
                rc = (rc << 2) | uint64_t((tmp & 3) ^ 3);
        }

        iceberg_upsert_increment(&table, fwd < rc ? fwd : rc, 0);
        ++inserted;

        // ── Slide remaining k-mers ───────────────────────────────────────────
        const uint8_t* byte_ptr = packed + (k >> 2);
        int            shift    = static_cast<int>(6 - 2 * (k & 3u));

        for (size_t i = k; i < len; ++i) {
            const uint64_t b = (*byte_ptr >> shift) & 3u;
            shift -= 2;
            if (shift < 0) { shift = 6; ++byte_ptr; }

            fwd = ((fwd << 2) | b) & mask;
            rc  = (rc >> 2) | ((b ^ 3) << (2 * (k - 1)));

            iceberg_upsert_increment(&table, fwd < rc ? fwd : rc, 0);
            ++inserted;
        }
    }

    if (dbg) {
        dbg->n_inserted = inserted;
        dbg->n_unique   = tot_balls(&table);
        dbg->table_cap  = table.metadata.nslots;
    }

    return inserted;
}


// ─── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct WriteCtx {
    const Config&  cfg;
    std::string&   chunk;
    std::ofstream& out;
    std::mutex&    out_mutex;
    uint64_t       written = 0;
};

template <uint16_t k, uint16_t m>
static void write_cb(KeyType key, ValueType cnt, void *ctx_ptr) {
    constexpr size_t WRITE_BATCH = 1u << 20;
    auto& ctx = *static_cast<WriteCtx<k, m>*>(ctx_ptr);
    if (cnt < ctx.cfg.ci || cnt > ctx.cfg.cx) return;

    std::string& chunk = ctx.chunk;
    std::string  label;
    kmer_to_str<k>(key, label);
    chunk += label;
    chunk += '\t';
    char cnt_buf[32];
    const auto [ptr, ec] = std::to_chars(std::begin(cnt_buf), std::end(cnt_buf),
                                          static_cast<uint32_t>(cnt));
    if (ec == std::errc()) chunk.append(cnt_buf, ptr);
    else chunk += std::to_string(cnt);
    chunk += '\n';
    ++ctx.written;

    if (chunk.size() >= WRITE_BATCH) {
        std::lock_guard<std::mutex> g(ctx.out_mutex);
        ctx.out.write(chunk.data(), chunk.size());
        chunk.clear();
    }
}

template <uint16_t k, uint16_t m>
uint64_t write_counts(
    iceberg_table& table,
    const Config&  cfg,
    std::string&   chunk,
    std::ofstream& out,
    std::mutex&    out_mutex)
{
    if (chunk.capacity() < (1u << 20) + 64) chunk.reserve((1u << 20) + 64);

    WriteCtx<k, m> ctx{cfg, chunk, out, out_mutex, 0};
    iceberg_iterate(&table,
                    static_cast<iceberg_iterate_cb>(write_cb<k, m>),
                    &ctx);
    if (!chunk.empty()) {
        std::lock_guard<std::mutex> g(out_mutex);
        out.write(chunk.data(), chunk.size());
        chunk.clear();
    }
    return ctx.written;
}


// ─── KFF output brick ─────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct KffCtx {
    const Config&         cfg;
    KffOutput&            kff_out;
    std::vector<uint8_t>& seq_buf;
    std::vector<uint32_t>& cnt_buf;
    uint64_t              written = 0;
};

template <uint16_t k, uint16_t m>
static void kff_cb(KeyType key, ValueType cnt, void *ctx_ptr) {
    constexpr size_t KMER_BYTES  = (k + 3) / 4;
    constexpr size_t BATCH_KMERS = (size_t(1) << 20) / (KMER_BYTES + 4);

    auto& ctx = *static_cast<KffCtx<k, m>*>(ctx_ptr);
    if (cnt < ctx.cfg.ci || cnt > ctx.cfg.cx) return;

    ctx.seq_buf.resize(ctx.seq_buf.size() + KMER_BYTES, 0);
    kmer_to_packed_msb<k>(key, ctx.seq_buf.data() + ctx.seq_buf.size() - KMER_BYTES);
    ctx.cnt_buf.push_back(static_cast<uint32_t>(cnt));
    ++ctx.written;

    if (ctx.cnt_buf.size() >= BATCH_KMERS) {
        ctx.kff_out.write_batch(ctx.seq_buf.data(), ctx.cnt_buf.data(), ctx.cnt_buf.size());
        ctx.seq_buf.clear();
        ctx.cnt_buf.clear();
    }
}

template <uint16_t k, uint16_t m>
uint64_t write_counts_kff(
    iceberg_table& table,
    const Config&  cfg,
    KffOutput&     kff_out)
{
    constexpr size_t KMER_BYTES  = (k + 3) / 4;
    constexpr size_t BATCH_KMERS = (size_t(1) << 20) / (KMER_BYTES + 4);

    std::vector<uint8_t>  seq_buf;
    std::vector<uint32_t> cnt_buf;
    seq_buf.reserve(BATCH_KMERS * KMER_BYTES);
    cnt_buf.reserve(BATCH_KMERS);

    KffCtx<k, m> ctx{cfg, kff_out, seq_buf, cnt_buf, 0};
    iceberg_iterate(&table,
                    static_cast<iceberg_iterate_cb>(kff_cb<k, m>),
                    &ctx);
    if (!cnt_buf.empty())
        kff_out.write_batch(seq_buf.data(), cnt_buf.data(), cnt_buf.size());
    return ctx.written;
}


// ─── Callback output brick ────────────────────────────────────────────────────

template <uint16_t k, uint16_t m, typename Callback>
struct CbCtx {
    const Config& cfg;
    Callback&     cb;
    std::string   label;
    uint64_t      written = 0;
};

template <uint16_t k, uint16_t m, typename Callback>
static void user_cb(KeyType key, ValueType cnt, void *ctx_ptr) {
    auto& ctx = *static_cast<CbCtx<k, m, Callback>*>(ctx_ptr);
    if (cnt < ctx.cfg.ci || cnt > ctx.cfg.cx) return;
    kmer_to_str<k>(key, ctx.label);
    ctx.cb(std::string_view(ctx.label), static_cast<uint32_t>(cnt));
    ++ctx.written;
}

template <uint16_t k, uint16_t m, typename Callback>
uint64_t write_counts_callback(
    iceberg_table& table,
    const Config&  cfg,
    Callback&      cb)
{
    CbCtx<k, m, Callback> ctx{cfg, cb, {}, 0};
    iceberg_iterate(&table,
                    [](KeyType key, ValueType cnt, void *p) {
                        user_cb<k, m, Callback>(key, cnt, p);
                    },
                    &ctx);
    return ctx.written;
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
    uint64_t sum_unique = 0, min_unique = UINT64_MAX, max_unique = 0;
    uint64_t sum_init   = 0, min_init   = UINT64_MAX, max_init   = 0;
    size_t   n_valid = 0;

    for (size_t p = 0; p < n_parts; ++p) {
        const auto& d = part_infos[p];
        if (d.table_cap == 0) continue;
        ++n_valid;
        const double lf = static_cast<double>(d.n_unique) / d.table_cap;
        sum_lf     += lf;
        min_lf      = std::min(min_lf, lf);
        max_lf      = std::max(max_lf, lf);
        sum_init   += d.init_sz;
        min_init    = std::min(min_init, d.init_sz);
        max_init    = std::max(max_init, d.init_sz);
        sum_unique += d.n_unique;
        min_unique  = std::min(min_unique, d.n_unique);
        max_unique  = std::max(max_unique, d.n_unique);
    }

    if (n_valid > 0) {
        const double mean_lf     = sum_lf     / n_valid;
        const double mean_unique = static_cast<double>(sum_unique) / n_valid;
        const double mean_init   = static_cast<double>(sum_init)   / n_valid;
        std::cerr << "\n[debug] aggregate table stats (" << n_parts << " partitions):\n"
                  << std::fixed;
        std::cerr << "  load_factor:   mean=" << std::setprecision(3) << mean_lf
                  << "  min=" << min_lf << "  max=" << max_lf << "\n";
        std::cerr << "  n_unique/part: mean=" << std::setprecision(0) << mean_unique
                  << "  min=" << min_unique << "  max=" << max_unique << "\n";
        std::cerr << "  init_sz/part:  mean=" << mean_init
                  << "  min=" << min_init   << "  max=" << max_init   << "\n";
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
            csv << p           << ","
                << d.init_sz   << ","
                << d.table_cap << ","
                << d.n_inserted << ","
                << d.n_unique   << ","
                << std::fixed << std::setprecision(6) << lf << "\n";
        }
        std::cerr << "[debug] table stats CSV: " << csv_path << "\n";
    }
}


// ─── Parallel harness (disk) ──────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    uint64_t       total_kmers,
    std::ofstream* out,
    KffOutput*     kff_out)
{
    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t> next_part{0};

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);

    auto worker = [&](size_t /*tid*/) {
        try {
            std::string chunk;

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;

                const std::string path = partition_path(cfg.work_dir, p);
                SuperkmerReader<k, m> reader(path);
                if (!reader.ok())
                    throw std::runtime_error("tuna: cannot open partition file for reading: " + path);

                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal)),
                                         size_t(1u << 15), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }

                iceberg_table table;
                iceberg_init(&table, iceberg_log_slots(init_sz));

                PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;
                if (dbg) { dbg->init_sz = init_sz; dbg->table_cap = table.metadata.nslots; }

                const uint64_t ins = count_partition<k, m>(reader, table, dbg);
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                if (cal == 0) {
                    const uint64_t unique = tot_balls(&table);
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                const uint64_t wrt = kff_out
                    ? write_counts_kff<k, m>(table, cfg, *kff_out)
                    : write_counts<k, m>(table, cfg, chunk, *out, out_mutex);
                total_written.fetch_add(wrt, std::memory_order_relaxed);

                iceberg_destroy(&table);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);

    if (cfg.debug_stats)
        emit_debug_stats(part_infos, n_parts, cfg);

    return { total_inserted.load(), total_written.load() };
}


// ─── In-memory counting harness ───────────────────────────────────────────────

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
    std::atomic<size_t> next_part{0};

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};

    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        std::string chunk;

        while (true) {
            const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
            if (p >= n_parts) break;

            const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
            size_t init_sz;
            if (cal > 0) {
                init_sz = std::clamp(next_pow2(static_cast<size_t>(cal)),
                                     size_t(1u << 15), size_t(1u << 22));
            } else {
                const size_t per_part = (total_kmers > 0)
                    ? static_cast<size_t>(total_kmers / n_parts) * 2
                    : size_t(1u << 27) / n_parts;
                init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
            }

            iceberg_table table;
            iceberg_init(&table, iceberg_log_slots(init_sz));

            PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;
            if (dbg) { dbg->init_sz = init_sz; dbg->table_cap = table.metadata.nslots; }

            uint64_t ins;
            {
                MemoryReader<k, m> reader(part_bufs[p]);
                ins = count_partition<k, m, MemoryReader<k, m>>(reader, table, dbg);
            }
            { std::string tmp; part_bufs[p].swap(tmp); }
            total_inserted.fetch_add(ins, std::memory_order_relaxed);

            if (cal == 0) {
                const uint64_t unique = tot_balls(&table);
                if (unique > 0) {
                    uint64_t expected = 0;
                    calibrated_unique.compare_exchange_strong(
                        expected, unique,
                        std::memory_order_relaxed, std::memory_order_relaxed);
                }
            }

            const uint64_t wrt = kff_out
                ? write_counts_kff<k, m>(table, cfg, *kff_out)
                : write_counts<k, m>(table, cfg, chunk, *out, out_mutex);
            total_written.fetch_add(wrt, std::memory_order_relaxed);

            iceberg_destroy(&table);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    if (cfg.debug_stats)
        emit_debug_stats(part_infos, n_parts, cfg);

    return { total_inserted.load(), total_written.load() };
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
    std::atomic<size_t> next_part{0};

    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        try {
            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;

                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal)),
                                         size_t(1u << 15), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }

                iceberg_table table;
                iceberg_init(&table, iceberg_log_slots(init_sz));

                uint64_t ins;
                {
                    MemoryReader<k, m> reader(part_bufs[p]);
                    ins = count_partition<k, m, MemoryReader<k, m>>(reader, table);
                }
                { std::string tmp; part_bufs[p].swap(tmp); }
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                if (cal == 0) {
                    const uint64_t unique = tot_balls(&table);
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);
                total_written.fetch_add(wrt, std::memory_order_relaxed);

                iceberg_destroy(&table);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
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
    std::atomic<size_t> next_part{0};

    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool>     stop{false};
    std::exception_ptr    worker_error = nullptr;
    std::mutex            worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        try {
            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;

                const std::string path = partition_path(cfg.work_dir, p);
                SuperkmerReader<k, m> reader(path);
                if (!reader.ok())
                    throw std::runtime_error("tuna: cannot open partition file for reading: " + path);

                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal)),
                                         size_t(1u << 15), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }

                iceberg_table table;
                iceberg_init(&table, iceberg_log_slots(init_sz));

                const uint64_t ins = count_partition<k, m>(reader, table);
                total_inserted.fetch_add(ins, std::memory_order_relaxed);

                if (cal == 0) {
                    const uint64_t unique = tot_balls(&table);
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);
                total_written.fetch_add(wrt, std::memory_order_relaxed);

                iceberg_destroy(&table);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
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
