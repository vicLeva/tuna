#pragma once

// Phase 2+3 — partition files → counted k-mers → output.
//
// Three independently modifiable bricks:
//
//   count_partition<k,m>   — read one partition file, fill a hash table.
//   write_counts<k,m>      — drain one table to the output stream.
//   count_and_write<k,m>   — parallel harness (threading / scheduling).

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
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <exception>
#include <stdexcept>
#include <charconv>
#include <cstring>
#include <string_view>
#include <filesystem>

#include <ankerl/unordered_dense.h>
#include "xxHash/xxhash.h"
#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"
#include "dense_kmer_table.hpp"

// Phase-2 counting table. Build with -DTUNA_PHASE2_DENSE to swap the kache
// quotienting table for an ankerl::unordered_dense map; see
// include/dense_kmer_table.hpp for what that costs and what it buys.
#ifdef TUNA_PHASE2_DENSE
template <uint16_t k_, bool mt_, typename T_, uint16_t l_, bool canonical_>
using phase2_table_t = tuna::Dense_Kmer_Hash_Table<k_, mt_, T_, l_, canonical_>;
#else
template <uint16_t k_, bool mt_, typename T_, uint16_t l_, bool canonical_>
using phase2_table_t = kache_hash::Streaming_Kmer_Hash_Table<k_, mt_, T_, l_, canonical_>;
#endif


// Returns the smallest power of two >= v (returns 1 for v == 0).
inline size_t next_pow2(size_t v) noexcept {
    size_t p = 1;
    while (p < v) p <<= 1;
    return p;
}


// ─── Debug stats ──────────────────────────────────────────────────────────────
// Populated by count_partition when dbg != nullptr (-dbg flag).
// coverage_hist is aggregated globally and written to debug_min_coverage.csv.

struct PartitionDebugInfo {
    uint64_t init_sz     = 0;   // initial hash table capacity (slots) at construction
    uint64_t n_inserted  = 0;   // total k-mer insertions (with multiplicity)
    uint64_t n_unique    = 0;   // unique k-mers stored (table.size() after counting)
    uint64_t n_overflow  = 0;   // insertions that went to overflow table
    uint64_t table_cap   = 0;   // final flat-table capacity (k-mer slots, after any resizes)
    uint64_t n_buckets   = 0;   // table_cap / B

    // Resize events recorded during count_partition for this partition.
    std::vector<kache_hash::ResizeEvent> resize_log;

    // coverage → count_of_minimizers_with_that_coverage
    // Only populated for superkmers with min_pos != 0xFFFF.
    std::unordered_map<uint32_t, uint64_t> coverage_hist;
};


// ─── Counting brick ───────────────────────────────────────────────────────────
//
// Drains a superkmer reader into the caller-provided hash table (mt_=false).
// Two windows alternate so the next superkmer is initialized and prefetched
// while the current one is still available to hide the bucket lookup latency.
//
// Returns the total number of k-mer insertions (with multiplicity).

template <uint16_t k, uint16_t m, bool canonical_ = true, typename Reader = SuperkmerReader<k, m>>
uint64_t count_partition(
    Reader&                                                                          reader,
    phase2_table_t<k, false, uint32_t, m, canonical_>&        table,
    typename phase2_table_t<k, false, uint32_t, m, canonical_>::Token& token,
    PartitionDebugInfo* dbg = nullptr)
{
    using hdr_t = sk_hdr_t<k, m>;              // superkmer header type (local alias)
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;  // sentinel: no precomputed minimizer hash

    auto inc = [](uint32_t v) { return v + 1; };
    uint64_t inserted = 0;

    // Per-minimizer k-mer count (only allocated when debug is requested).
    std::unordered_map<uint64_t, uint32_t> min_kmer_count;

    // ── Prime the pump ────────────────────────────────────────────────────────
    if (!reader.next()) return inserted;

    const uint8_t* cur_packed  = reader.packed_data();  // packed bases of current superkmer
    size_t         cur_len     = reader.size();
    hdr_t          cur_min_pos = reader.min_pos();       // minimizer position in current superkmer

    std::array<kache_hash::Kmer_Window<k, m>, 2> windows;
    auto* win = &windows[0];
    auto* next_win = &windows[1];
    if (cur_len >= k) {
        if (cur_min_pos != NO_MIN) {
            const uint64_t mh = win->init_packed_with_min(cur_packed, cur_min_pos);
            if (dbg) min_kmer_count[mh] += static_cast<uint32_t>(cur_len - k + 1);
        } else {
            win->init_packed_known_out(cur_packed);
        }
        table.prefetch(*win);
    }

    // ── Main loop ─────────────────────────────────────────────────────────────
    while (reader.next()) {
        const uint8_t* nxt_packed  = reader.packed_data();  // packed bases of next superkmer (prefetch target)
        const size_t   nxt_len     = reader.size();
        const hdr_t    nxt_min_pos = reader.min_pos();

        if (nxt_len >= k) {
            if (nxt_min_pos != NO_MIN) {
                const uint64_t mh = next_win->init_packed_with_min(nxt_packed, nxt_min_pos);
                if (dbg) min_kmer_count[mh] += static_cast<uint32_t>(nxt_len - k + 1);
            } else {
                next_win->init_packed_known_out(nxt_packed);
            }
            table.prefetch(*next_win);
        }

        // Process CURRENT superkmer.
        if (cur_len >= k) {
            table.upsert(*win, inc, uint32_t(1), token);
            ++inserted;

            // Unpack subsequent bases directly as DNA::Base (kache encoding).
            const uint8_t* byte_ptr = cur_packed + (k >> 2);
            int shift = static_cast<int>(6u - 2u * (k & 3u));
            for (size_t i = k; i < cur_len; ++i) {
                const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
                shift -= 2;
                if (shift < 0) { shift = 6; ++byte_ptr; }
                win->advance_known_out(b);
                table.upsert(*win, inc, uint32_t(1), token);
                ++inserted;
            }
        }

        // Advance to the initialized and already-prefetched window.
        cur_packed  = nxt_packed;
        cur_len     = nxt_len;
        std::swap(win, next_win);
    }

    // ── Last superkmer ────────────────────────────────────────────────────────
    if (cur_len >= k) {
        table.upsert(*win, inc, uint32_t(1), token);
        ++inserted;
        const uint8_t* byte_ptr = cur_packed + (k >> 2);
        int shift = static_cast<int>(6u - 2u * (k & 3u));
        for (size_t i = k; i < cur_len; ++i) {
            const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
            shift -= 2;
            if (shift < 0) { shift = 6; ++byte_ptr; }
            win->advance_known_out(b);
            table.upsert(*win, inc, uint32_t(1), token);
            ++inserted;
        }
    }

    if (dbg) {
        // Convert per-minimizer k-mer counts to a coverage histogram.
        for (auto& [mh, cnt] : min_kmer_count)
            dbg->coverage_hist[cnt]++;
    }

    return inserted;
}

// Insert every k-mer represented by one packed superkmer, increasing counts by
// multiplicity. This is the replay primitive for exact record aggregation.
template <uint16_t k, uint16_t m, bool canonical_ = true>
uint64_t count_initialized_superkmer_record(
    const uint8_t* packed,
    size_t len,
    uint32_t multiplicity,
    kache_hash::Kmer_Window<k, m>& win,
    phase2_table_t<k, false, uint32_t, m, canonical_>& table,
    typename phase2_table_t<k, false, uint32_t, m, canonical_>::Token& token)
{
    if (len < k || multiplicity == 0) return 0;

    auto inc = [multiplicity](uint32_t v) { return v + multiplicity; };
    table.upsert(win, inc, multiplicity, token);

    const uint8_t* byte_ptr = packed + (k >> 2);
    int shift = static_cast<int>(6u - 2u * (k & 3u));
    for (size_t i = k; i < len; ++i) {
        const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
        shift -= 2;
        if (shift < 0) { shift = 6; ++byte_ptr; }
        win.advance_known_out(b);
        table.upsert(win, inc, multiplicity, token);
    }

    return static_cast<uint64_t>(multiplicity) * (len - k + 1);
}

template <uint16_t k, uint16_t m, bool canonical_ = true>
uint64_t count_superkmer_record(
    const uint8_t* packed,
    size_t len,
    uint32_t multiplicity,
    phase2_table_t<k, false, uint32_t, m, canonical_>& table,
    typename phase2_table_t<k, false, uint32_t, m, canonical_>::Token& token)
{
    if (len < k || multiplicity == 0) return 0;

    kache_hash::Kmer_Window<k, m> win;
    win.init_packed_known_out(packed);
    table.prefetch(win);
    return count_initialized_superkmer_record<k, m, canonical_>(
        packed, len, multiplicity, win, table, token);
}

inline size_t superkmer_dedup_worker_budget(const Config& cfg, size_t n_threads)
{
    static constexpr uint64_t MiB = uint64_t(1) << 20;
    if (cfg.ram_budget_bytes == 0) return static_cast<size_t>(256 * MiB);

    const uint64_t share =
        (cfg.ram_budget_bytes * 3u / 8u) / std::max<uint64_t>(1, n_threads);
    return static_cast<size_t>(
        std::clamp<uint64_t>(share, 16 * MiB, 512 * MiB));
}

struct SuperkmerRecordHash {
    using is_avalanching = void;

    size_t operator()(std::string_view record) const noexcept
    {
        return static_cast<size_t>(XXH3_64bits(record.data(), record.size()));
    }
};

struct SuperkmerDedupStats {
    std::atomic<uint64_t> records{0};
    std::atomic<uint64_t> unique_records{0};
    std::atomic<uint64_t> aggregate_ns{0};
    std::atomic<uint64_t> replay_ns{0};
    std::atomic<uint64_t> aggregated_partitions{0};
    std::atomic<uint64_t> direct_partitions{0};
};

// Aggregate identical encoded records while their mmap/string backing remains
// alive, then replay each record once with its exact multiplicity. In auto mode
// a low-redundancy sample or a memory-limit hit resets the reader and uses the
// direct path before anything has been inserted into the k-mer table.
template <uint16_t k, uint16_t m, bool canonical_ = true, typename Reader>
uint64_t count_partition_aggregated(
    Reader& reader,
    SuperkmerDedupMode mode,
    size_t memory_budget,
    phase2_table_t<k, false, uint32_t, m, canonical_>& table,
    typename phase2_table_t<k, false, uint32_t, m, canonical_>::Token& token,
    PartitionDebugInfo* dbg = nullptr,
    SuperkmerDedupStats* stats = nullptr)
{
    using map_t =
        ankerl::unordered_dense::map<std::string_view, uint32_t, SuperkmerRecordHash>;
    static constexpr size_t ENTRY_ESTIMATE = 32;
    static constexpr size_t AUTO_SAMPLE_RECORDS = 16u << 10;
    static constexpr size_t AUTO_MIN_RECORDS = 4u << 10;
    static constexpr size_t AUTO_DUPLICATE_DENOMINATOR = 16;

    map_t multiplicities;
    const size_t reserve_by_data =
        std::max<size_t>(1024, reader.data_size() / 8);
    const size_t reserve_by_budget =
        std::max<size_t>(1024, memory_budget / ENTRY_ESTIMATE);
    multiplicities.reserve(std::min(reserve_by_data, reserve_by_budget));

    size_t records = 0;
    bool use_aggregation = true;
    const auto aggregate_start = std::chrono::steady_clock::now();
    while (reader.next()) {
        const std::string_view record(reader.record_data(), reader.record_size());
        auto [it, inserted] = multiplicities.try_emplace(record, uint32_t(1));
        if (!inserted) ++it->second;
        ++records;

        if ((records & 4095u) == 0) {
            if (multiplicities.size() * ENTRY_ESTIMATE > memory_budget) {
                use_aggregation = false;
                break;
            }
            if (mode == SuperkmerDedupMode::Auto &&
                records == AUTO_SAMPLE_RECORDS) {
                const size_t duplicates = records - multiplicities.size();
                if (duplicates * AUTO_DUPLICATE_DENOMINATOR < records) {
                    use_aggregation = false;
                    break;
                }
            }
        }
    }

    if (mode == SuperkmerDedupMode::Auto && records < AUTO_SAMPLE_RECORDS) {
        const size_t duplicates = records - multiplicities.size();
        use_aggregation =
            records >= AUTO_MIN_RECORDS &&
            duplicates * AUTO_DUPLICATE_DENOMINATOR >= records;
    }
    const auto aggregate_ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - aggregate_start).count());

    if (!use_aggregation) {
        if (stats) {
            stats->aggregate_ns.fetch_add(aggregate_ns, std::memory_order_relaxed);
            stats->direct_partitions.fetch_add(1, std::memory_order_relaxed);
        }
        multiplicities.clear();
        multiplicities.rehash(0);
        reader.reset();
        return count_partition<k, m, canonical_, Reader>(reader, table, token, dbg);
    }

    using hdr_t = sk_hdr_t<k, m>;
    uint64_t total_inserted = 0;
    const auto replay_start = std::chrono::steady_clock::now();

    auto decode_record = [](const std::string_view record,
                            const uint8_t*& packed,
                            size_t& len) {
        hdr_t encoded_len;
        std::memcpy(&encoded_len, record.data(), sizeof(hdr_t));
        len = static_cast<size_t>(encoded_len);
        packed = reinterpret_cast<const uint8_t*>(
            record.data() + sizeof(hdr_t));
    };

    auto it = multiplicities.begin();
    const auto end = multiplicities.end();
    std::array<kache_hash::Kmer_Window<k, m>, 2> windows;
    auto* win = &windows[0];
    auto* next_win = &windows[1];

    const uint8_t* cur_packed = nullptr;
    size_t cur_len = 0;
    uint32_t cur_multiplicity = 0;
    if (it != end) {
        decode_record(it->first, cur_packed, cur_len);
        cur_multiplicity = it->second;
        win->init_packed_known_out(cur_packed);
        table.prefetch(*win);
        ++it;
    }

    while (it != end) {
        const uint8_t* next_packed = nullptr;
        size_t next_len = 0;
        decode_record(it->first, next_packed, next_len);
        const uint32_t next_multiplicity = it->second;
        next_win->init_packed_known_out(next_packed);
        table.prefetch(*next_win);

        total_inserted += count_initialized_superkmer_record<k, m, canonical_>(
            cur_packed, cur_len, cur_multiplicity, *win, table, token);

        cur_packed = next_packed;
        cur_len = next_len;
        cur_multiplicity = next_multiplicity;
        std::swap(win, next_win);
        ++it;
    }

    if (cur_packed != nullptr) {
        total_inserted += count_initialized_superkmer_record<k, m, canonical_>(
            cur_packed, cur_len, cur_multiplicity, *win, table, token);
    }
    if (stats) {
        const auto replay_ns = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - replay_start).count());
        stats->records.fetch_add(records, std::memory_order_relaxed);
        stats->unique_records.fetch_add(multiplicities.size(), std::memory_order_relaxed);
        stats->aggregate_ns.fetch_add(aggregate_ns, std::memory_order_relaxed);
        stats->replay_ns.fetch_add(replay_ns, std::memory_order_relaxed);
        stats->aggregated_partitions.fetch_add(1, std::memory_order_relaxed);
    }
    return total_inserted;
}

inline void emit_superkmer_dedup_stats(const SuperkmerDedupStats& stats)
{
    const uint64_t records = stats.records.load(std::memory_order_relaxed);
    const uint64_t unique = stats.unique_records.load(std::memory_order_relaxed);
    const double multiplicity =
        unique > 0 ? static_cast<double>(records) / static_cast<double>(unique) : 0.0;
    std::cerr << "[dedup] records=" << records
              << " unique_records=" << unique
              << " multiplicity=" << std::fixed << std::setprecision(3) << multiplicity
              << " aggregate_cpu_s="
              << (static_cast<double>(stats.aggregate_ns.load(std::memory_order_relaxed)) / 1e9)
              << " replay_cpu_s="
              << (static_cast<double>(stats.replay_ns.load(std::memory_order_relaxed)) / 1e9)
              << " aggregated_partitions="
              << stats.aggregated_partitions.load(std::memory_order_relaxed)
              << " direct_partitions="
              << stats.direct_partitions.load(std::memory_order_relaxed) << "\n";
}


// ─── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m, bool mt_ = false, bool canonical_ = true>
uint64_t write_counts(
    phase2_table_t<k, mt_, uint32_t, m, canonical_>& table,
    const Config&   cfg,
    std::string&    chunk,
    std::ofstream&  out,
    std::mutex&     out_mutex)
{
    constexpr size_t WRITE_BATCH = 1u << 20; // 1 MB
    if (chunk.capacity() < WRITE_BATCH + 64) chunk.reserve(WRITE_BATCH + 64);

    const auto flush_chunk = [&]() {
        std::lock_guard<std::mutex> g(out_mutex);
        out.write(chunk.data(), chunk.size());
        chunk.clear();
    };

    std::string label;
    uint64_t written = 0;

    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        if (cnt < cfg.ci || cnt > cfg.cx) return;
        entry.first.get_label(label);
        chunk += label;
        chunk += '\t';
        char cnt_buf[32];
        const auto [ptr, ec] = std::to_chars(std::begin(cnt_buf), std::end(cnt_buf), cnt);
        if (ec == std::errc()) chunk.append(cnt_buf, ptr);
        else chunk += std::to_string(cnt);
        chunk += '\n';
        ++written;
        if (chunk.size() >= WRITE_BATCH) flush_chunk();
    });
    if (!chunk.empty()) flush_chunk();
    return written;
}


// ─── KFF output brick ─────────────────────────────────────────────────────────
//
// Encodes complete KFF records in four worker-local buffers, one per lossless
// count width. Batching by width avoids a compaction pass and keeps section
// changes bounded while giving every record its minimum 1-4 byte count.

template <uint16_t k>
class KffBatchWriter
{
public:
    static constexpr size_t KMER_BYTES  = (k + 3) / 4;
    static constexpr size_t BATCH_BYTES = size_t(1) << 20;

    explicit KffBatchWriter(KffOutput* out)
        : out_(out)
    {
        if (out_) {
            init_buffer<1>();
            init_buffer<2>();
            init_buffer<3>();
            init_buffer<4>();
        }
    }

    template <typename Table>
    uint64_t write(Table& table, const Config& cfg)
    {
        uint64_t written = 0;
        table.for_each([&](const auto& entry) {
            const uint64_t cnt = entry.second;
            if (cnt < cfg.ci || cnt > cfg.cx) return;

            const uint32_t c = static_cast<uint32_t>(cnt);
            if (c <= std::numeric_limits<uint8_t>::max())
                append<1>(entry.first, c);
            else if (c <= std::numeric_limits<uint16_t>::max())
                append<2>(entry.first, c);
            else if (c <= 0xFFFFFFu)
                append<3>(entry.first, c);
            else
                append<4>(entry.first, c);
            ++written;
        });
        return written;
    }

    void flush()
    {
        flush_buffer<1>();
        flush_buffer<2>();
        flush_buffer<3>();
        flush_buffer<4>();
    }

private:
    struct Buffer {
        std::vector<uint8_t> records;
        size_t buffered = 0;
    };

    template <size_t COUNT_BYTES>
    static constexpr size_t batch_kmers()
    {
        return BATCH_BYTES / (KMER_BYTES + COUNT_BYTES);
    }

    template <size_t COUNT_BYTES>
    void init_buffer()
    {
        auto& buffer = buffers_[COUNT_BYTES - 1];
        buffer.records.resize(
            batch_kmers<COUNT_BYTES>() * (KMER_BYTES + COUNT_BYTES));
    }

    template <size_t COUNT_BYTES, typename Kmer>
    void append(const Kmer& kmer, uint32_t count)
    {
        constexpr size_t RECORD_SIZE = KMER_BYTES + COUNT_BYTES;
        auto& buffer = buffers_[COUNT_BYTES - 1];
        uint8_t* dst = buffer.records.data() + buffer.buffered * RECORD_SIZE;
        kmer.write_packed_2bit_msb(dst);
        for (size_t i = 0; i < COUNT_BYTES; ++i) {
            const size_t shift = 8 * (COUNT_BYTES - 1 - i);
            dst[KMER_BYTES + i] = static_cast<uint8_t>(count >> shift);
        }
        if (++buffer.buffered >= batch_kmers<COUNT_BYTES>())
            flush_buffer<COUNT_BYTES>();
    }

    template <size_t COUNT_BYTES>
    void flush_buffer()
    {
        auto& buffer = buffers_[COUNT_BYTES - 1];
        if (buffer.buffered == 0) return;
        out_->write_batch(
            buffer.records.data(), buffer.buffered, COUNT_BYTES);
        buffer.buffered = 0;
    }

    KffOutput* out_;
    std::array<Buffer, 4> buffers_;
};

template <uint16_t k, uint16_t m, bool mt_ = false, bool canonical_ = true>
uint64_t count_filtered(
    const phase2_table_t<k, mt_, uint32_t, m, canonical_>& table,
    const Config& cfg)
{
    if (cfg.ci <= 1 && cfg.cx >= std::numeric_limits<uint32_t>::max())
        return static_cast<uint64_t>(table.size());

    uint64_t written = 0;
    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        written += cnt >= cfg.ci && cnt <= cfg.cx;
    });
    return written;
}


// ─── Callback output brick ────────────────────────────────────────────────────
//
// Drains the table, applies ci/cx filters, and calls cb(kmer, count) for each
// passing k-mer.  cb is called from the calling thread only (one table per call),
// so no internal mutex is used.  If multiple worker threads drain different tables
// concurrently, cb may be invoked from several threads simultaneously — the caller
// is responsible for any needed synchronisation.

template <uint16_t k, uint16_t m, bool mt_ = false, bool canonical_ = true, typename Callback>
uint64_t write_counts_callback(
    phase2_table_t<k, mt_, uint32_t, m, canonical_>& table,
    const Config& cfg,
    Callback& cb)
{
    uint64_t written = 0;
    [[maybe_unused]] std::string label;

    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        if (cnt < cfg.ci || cnt > cfg.cx) return;
        if constexpr (std::is_invocable_v<Callback, const kache_hash::Kmer<k>&, uint32_t>) {
            cb(entry.first, static_cast<uint32_t>(cnt));
        } else {
            entry.first.get_label(label);
            cb(std::string_view(label), static_cast<uint32_t>(cnt));
        }
        ++written;
    });

    return written;
}


// ─── Callback counting harnesses ──────────────────────────────────────────────
//
// Same structure as count_and_write / count_and_write_mem but drain tables via
// a user-supplied callback instead of writing to a file.
//
// Thread safety: cb may be called concurrently from multiple worker threads
// (one per partition).  Each partition's k-mers are disjoint, so there is no
// risk of duplicate calls for the same k-mer.  The caller must ensure cb is
// safe to call from multiple threads if num_threads > 1.

template <uint16_t k, uint16_t m, bool canonical_ = true, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    Callback&&                cb)
{
    using table_t = phase2_table_t<k, false, uint32_t, m, canonical_>;  // hash table type (local alias)

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t> next_part{0};

    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        try {
            typename table_t::Token token;

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                         size_t(1u << 16), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }
                table_t table(init_sz, 1);

                uint64_t ins;  // k-mers inserted into this partition
                {
                    MemoryReader<k, m> reader(part_bufs[p]);
                    ins = count_partition<k, m, canonical_, MemoryReader<k, m>>(reader, table, token);
                }
                { std::string tmp; part_bufs[p].swap(tmp); }
                total_inserted.fetch_add(ins, std::memory_order_relaxed);
                if (cal == 0) {
                    const uint64_t unique = static_cast<uint64_t>(table.size());
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);  // k-mers written to output
                total_written.fetch_add(wrt, std::memory_order_relaxed);
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


template <uint16_t k, uint16_t m, bool canonical_ = true, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback(
    const Config& cfg,
    uint64_t      total_kmers,
    Callback&&    cb)
{
    using table_t = phase2_table_t<k, false, uint32_t, m, canonical_>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t> next_part{0};

    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex worker_error_mutex;
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        try {
            typename table_t::Token token;

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const std::string path = partition_path(cfg.work_dir, p);
                SuperkmerReader<k, m> reader(path);
                if (!reader.ok())
                    throw std::runtime_error(
                        "tuna: cannot open partition file for reading: " + path);
                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                         size_t(1u << 16), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }
                table_t table(init_sz, 1);

                const uint64_t ins = count_partition<k, m, canonical_>(reader, table, token);
                total_inserted.fetch_add(ins, std::memory_order_relaxed);
                if (cal == 0) {
                    const uint64_t unique = static_cast<uint64_t>(table.size());
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                const uint64_t wrt = write_counts_callback<k, m>(table, cfg, cb);
                total_written.fetch_add(wrt, std::memory_order_relaxed);
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


// ─── Debug output helper ──────────────────────────────────────────────────────
//
// Called by count_and_write (disk) and count_and_write_mem (in-memory) after
// all partitions have been processed.  Prints per-partition table stats,
// aggregate load-factor analysis, resize summary, and minimizer coverage.
// Also writes to cfg.work_dir:
//   debug_table_stats.csv   — one row per partition (all partitions)
//   debug_min_coverage.csv  — minimizer coverage histogram

inline void emit_debug_stats(
    const std::vector<PartitionDebugInfo>& part_infos,
    size_t n_parts,
    const Config& cfg)
{
    // ── Per-partition table summary (first 20) ────────────────────────────────
    std::cerr << "\n[debug] per-partition table stats (first 20):\n";
    std::cerr << "  part   init_sz   table_cap   n_inserted    n_unique  load%  n_resizes\n";
    for (size_t p = 0; p < std::min(n_parts, size_t(20)); ++p) {
        const auto& d = part_infos[p];
        const double lf = d.table_cap ? 100.0 * d.n_unique / d.table_cap : 0.0;
        std::cerr << "  " << std::setw(5)  << p
                  << "  " << std::setw(8)  << d.init_sz
                  << "  " << std::setw(9)  << d.table_cap
                  << "  " << std::setw(11) << d.n_inserted
                  << "  " << std::setw(10) << d.n_unique
                  << "  " << std::fixed << std::setprecision(1) << std::setw(5) << lf << "%"
                  << "  " << d.resize_log.size() << "\n";
    }
    if (n_parts > 20) std::cerr << "  ... (" << (n_parts - 20) << " more partitions)\n";

    // ── Aggregate stats ───────────────────────────────────────────────────────
    double   sum_lf = 0.0,                    min_lf = 2.0,                    max_lf = -1.0;
    double   sum_or = 0.0,                    min_or = 1e18,                   max_or = -1e18;
    uint64_t sum_unique = 0,                  min_unique = static_cast<uint64_t>(-1), max_unique = 0;
    uint64_t sum_init   = 0,                  min_init   = static_cast<uint64_t>(-1), max_init   = 0;
    size_t   n_valid = 0, n_nonempty = 0;
    uint64_t total_resizes = 0;
    size_t   n_parts_resized = 0;

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
        if (d.n_unique > 0) {
            ++n_nonempty;
            const double or_ = static_cast<double>(d.init_sz) / d.n_unique;
            sum_or += or_;
            min_or  = std::min(min_or, or_);
            max_or  = std::max(max_or, or_);
        }
        if (!d.resize_log.empty()) ++n_parts_resized;
        total_resizes += d.resize_log.size();
    }

    if (n_valid > 0) {
        const double mean_lf     = sum_lf / n_valid;
        const double mean_or     = n_nonempty ? sum_or    / n_nonempty : 0.0;
        const double mean_unique = static_cast<double>(sum_unique) / n_valid;
        const double mean_init   = static_cast<double>(sum_init)   / n_valid;

        std::cerr << "\n[debug] aggregate table stats (" << n_parts << " partitions):\n"
                  << std::fixed;
        std::cerr << "  load_factor:    mean=" << std::setprecision(3) << mean_lf
                  << "  min=" << min_lf << "  max=" << max_lf << "\n";
        if (n_nonempty > 0)
            std::cerr << "  oversize_ratio: mean=" << std::setprecision(1) << mean_or
                      << "  min=" << min_or << "  max=" << max_or
                      << "  (init_sz / n_unique)\n";
        std::cerr << "  n_unique/part:  mean=" << std::setprecision(0) << mean_unique
                  << "  min=" << min_unique << "  max=" << max_unique << "\n";
        std::cerr << "  init_sz/part:   mean=" << mean_init
                  << "  min=" << min_init << "  max=" << max_init << "\n";

        // Machine-parseable structured lines (grep-friendly for bench scripts).
        std::cerr << "dbg_load_mean: "     << std::setprecision(6) << mean_lf    << "\n"
                  << "dbg_load_min: "      << min_lf                              << "\n"
                  << "dbg_load_max: "      << max_lf                              << "\n"
                  << "dbg_oversize_mean: " << std::setprecision(2) << mean_or    << "\n"
                  << "dbg_unique_mean: "   << std::setprecision(0) << mean_unique << "\n"
                  << "dbg_n_resizes: "     << total_resizes                       << "\n"
                  << "dbg_parts_resized: " << n_parts_resized                     << "\n";
    }

    // ── Per-partition CSV (all partitions) ────────────────────────────────────
    {
        const std::string csv_path = cfg.work_dir + "debug_table_stats.csv";
        std::ofstream csv(csv_path);
        if (csv) {
            csv << "partition_id,init_sz,table_cap,n_inserted,n_unique,load_factor,n_resizes,resize_s\n";
            for (size_t p = 0; p < n_parts; ++p) {
                const auto& d = part_infos[p];
                const double lf = d.table_cap
                    ? static_cast<double>(d.n_unique) / d.table_cap : 0.0;
                double resize_s = 0.0;
                for (const auto& ev : d.resize_log) resize_s += ev.elapsed_s;
                csv << p             << ","
                    << d.init_sz     << ","
                    << d.table_cap   << ","
                    << d.n_inserted  << ","
                    << d.n_unique    << ","
                    << std::fixed << std::setprecision(6) << lf << ","
                    << d.resize_log.size() << ","
                    << std::setprecision(6) << resize_s << "\n";
            }
            std::cerr << "[debug] table stats CSV: " << csv_path << "\n";
        } else {
            std::cerr << "[debug] warning: could not write table stats CSV to " << csv_path << "\n";
        }
    }

    // ── Resize event summary ──────────────────────────────────────────────────
    {
        uint64_t total_resize_events = 0;
        double   total_resize_s      = 0.0;
        uint64_t n_ov_triggered      = 0;
        uint64_t n_load_triggered    = 0;
        size_t   n_parts_resized_    = 0;
        double   max_resize_s        = 0.0;
        size_t   max_resize_part     = 0;
        uint64_t max_resize_count    = 0;
        size_t   max_resize_count_part = 0;

        struct PartResizeSummary {
            size_t   part;
            uint64_t n_resizes;
            double   total_s;
            uint64_t n_ov;
            uint64_t n_load;
        };
        std::vector<PartResizeSummary> summaries;

        for (size_t p = 0; p < n_parts; ++p) {
            const auto& rlog = part_infos[p].resize_log;
            if (rlog.empty()) continue;
            ++n_parts_resized_;
            PartResizeSummary ps{p, rlog.size(), 0.0, 0, 0};
            for (const auto& ev : rlog) {
                ps.total_s += ev.elapsed_s;
                ps.n_ov    += ev.overflow_triggered ? 1 : 0;
                ps.n_load  += ev.overflow_triggered ? 0 : 1;
                if (ev.elapsed_s > max_resize_s) { max_resize_s = ev.elapsed_s; max_resize_part = p; }
            }
            total_resize_events += ps.n_resizes;
            total_resize_s      += ps.total_s;
            n_ov_triggered      += ps.n_ov;
            n_load_triggered    += ps.n_load;
            if (ps.n_resizes > max_resize_count) { max_resize_count = ps.n_resizes; max_resize_count_part = p; }
            summaries.push_back(ps);
        }

        std::cerr << "\n[debug] resize summary:\n";
        std::cerr << "  partitions with resizes : " << n_parts_resized_ << " / " << n_parts << "\n";
        std::cerr << "  total resize events     : " << total_resize_events << "\n";
        std::cerr << "  total resize time       : " << std::fixed << std::setprecision(3) << total_resize_s << "s\n";
        std::cerr << "  overflow-triggered      : " << n_ov_triggered << "\n";
        std::cerr << "  load-triggered          : " << n_load_triggered << "\n";
        if (total_resize_events > 0) {
            std::cerr << "  slowest single resize   : " << std::setprecision(3) << max_resize_s
                      << "s (part " << max_resize_part << ")\n";
            std::cerr << "  most resizes in one part: " << max_resize_count
                      << " (part " << max_resize_count_part << ")\n";
        }

        if (!summaries.empty()) {
            std::partial_sort(summaries.begin(),
                              summaries.begin() + std::min(summaries.size(), size_t(10)),
                              summaries.end(),
                              [](const auto& a, const auto& b){ return a.total_s > b.total_s; });
            const size_t show = std::min(summaries.size(), size_t(10));
            std::cerr << "\n  top " << show << " partitions by resize cost:\n";
            std::cerr << "  part   n_resizes  resize_s   n_ov_trig  n_load_trig\n";
            for (size_t i = 0; i < show; ++i) {
                const auto& ps = summaries[i];
                std::cerr << "  " << std::setw(5) << ps.part
                          << "  " << std::setw(9) << ps.n_resizes
                          << "  " << std::setprecision(3) << std::setw(9) << ps.total_s
                          << "  " << std::setw(9) << ps.n_ov
                          << "  " << std::setw(11) << ps.n_load << "\n";
                for (const auto& ev : part_infos[ps.part].resize_log)
                    std::cerr << "           "
                              << ev.old_cap << "->" << (ev.old_cap * 2)
                              << "  " << (ev.overflow_triggered ? "ov  " : "load")
                              << "  " << std::setprecision(3) << ev.elapsed_s << "s"
                              << "  ov_count=" << ev.ov_count
                              << "  main_sz=" << ev.main_sz << "\n";
            }
        }
    }

    // ── Minimizer coverage histogram summary ──────────────────────────────────
    {
        // Aggregate from per-partition histograms (populated by count_partition when dbg != nullptr).
        std::unordered_map<uint32_t, uint64_t> global_coverage_hist;
        for (size_t p = 0; p < n_parts; ++p)
            for (const auto& [cov, cnt] : part_infos[p].coverage_hist)
                global_coverage_hist[cov] += cnt;

        if (!global_coverage_hist.empty()) {
            uint64_t total_minimizers = 0, total_kmers_covered = 0;
            uint32_t max_cov = 0;
            for (auto& [cov, cnt] : global_coverage_hist) {
                total_minimizers    += cnt;
                total_kmers_covered += static_cast<uint64_t>(cov) * cnt;
                if (cov > max_cov) max_cov = cov;
            }
            const double avg_cov = total_minimizers
                ? static_cast<double>(total_kmers_covered) / total_minimizers : 0.0;

            std::cerr << "\n[debug] minimizer coverage (k-mers sharing one minimizer):\n";
            std::cerr << "  unique minimizers tracked : " << total_minimizers << "\n";
            std::cerr << "  total k-mers covered      : " << total_kmers_covered << "\n";
            std::cerr << "  avg k-mers / minimizer    : " << std::fixed << std::setprecision(2) << avg_cov << "\n";
            std::cerr << "  max k-mers / minimizer    : " << max_cov << "\n";

            const std::string csv_path = cfg.work_dir + "debug_min_coverage.csv";
            std::ofstream csv(csv_path);
            if (csv) {
                csv << "coverage,n_minimizers,total_kmers\n";
                std::vector<std::pair<uint32_t, uint64_t>> sorted(
                    global_coverage_hist.begin(), global_coverage_hist.end());
                std::sort(sorted.begin(), sorted.end(),
                          [](const auto& a, const auto& b){ return a.first < b.first; });
                for (auto& [cov, cnt] : sorted)
                    csv << cov << "," << cnt << ","
                        << (static_cast<uint64_t>(cov) * cnt) << "\n";
                std::cerr << "[debug] minimizer coverage CSV written to: " << csv_path << "\n";
            } else {
                std::cerr << "[debug] warning: could not write CSV to " << csv_path << "\n";
            }
        }
    }
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// Workers steal partitions round-robin (p = tid, tid+n_threads, …).
// init_sz = clamp(2 × total_kmers/n_parts, 4K, 4M).
// total_kmers counts all occurrences (with multiplicity), so 2× gives headroom
// for the unique fraction without oversizing beyond the 4M cap.

template <uint16_t k, uint16_t m, bool canonical_ = true>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    uint64_t       total_kmers,
    std::ofstream* out,       // non-null for TSV output
    KffOutput*     kff_out)   // non-null for KFF output
{
    using table_t = phase2_table_t<k, false, uint32_t, m, canonical_>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    const size_t dedup_budget = superkmer_dedup_worker_budget(cfg, n_threads);
    std::atomic<size_t> next_part{0};

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<bool> stop{false};
    std::exception_ptr worker_error = nullptr;
    std::mutex worker_error_mutex;

    // Overflow statistics accumulated across all partitions.
    std::mutex                                     ov_stats_mutex;
    std::vector<std::pair<uint32_t, uint64_t>>     ov_top_global;   // merged top minimizers
    std::atomic<uint64_t>                          ov_total{0};     // total k-mers sent to overflow table

    // Debug statistics (only used when cfg.debug_stats).
    std::vector<PartitionDebugInfo>                part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};
    SuperkmerDedupStats dedup_stats;

    auto worker = [&](size_t /*tid*/) {
        try {
            typename table_t::Token token;
            std::string chunk;
            KffBatchWriter<k> kff_writer(kff_out);

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const std::string path = partition_path(cfg.work_dir, p);
                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                         size_t(1u << 16), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }
                table_t table(init_sz, 1);

                PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;
                uint64_t ins;
                {
                    SuperkmerReader<k, m> reader(path);
                    if (!reader.ok())
                        throw std::runtime_error(
                            "tuna: cannot open partition file for reading: " + path);
                    ins =
                        cfg.dedup_mode == SuperkmerDedupMode::Off
                            ? count_partition<k, m, canonical_>(reader, table, token, dbg)
                            : count_partition_aggregated<k, m, canonical_>(
                                reader, cfg.dedup_mode, dedup_budget,
                                table, token, dbg,
                                cfg.debug_stats ? &dedup_stats : nullptr);
                }
                if (!cfg.keep_tmp) {
                    std::error_code ec;
                    std::filesystem::remove(path, ec);
                }
                total_inserted.fetch_add(ins, std::memory_order_relaxed);
                if (cal == 0) {
                    const uint64_t unique = static_cast<uint64_t>(table.size());
                    if (unique > 0) {
                        uint64_t expected = 0;
                        calibrated_unique.compare_exchange_strong(
                            expected, unique,
                            std::memory_order_relaxed, std::memory_order_relaxed);
                    }
                }

                if (dbg) {
                    dbg->init_sz     = init_sz;
                    dbg->n_inserted  = ins;
                    dbg->n_unique    = static_cast<uint64_t>(table.size());
                    dbg->n_overflow  = table.overflow_insert_count();
                    dbg->table_cap   = table.capacity();
                    dbg->n_buckets   = table.bucket_count();
                    dbg->resize_log  = table.resize_log();
                }

                const uint64_t wrt = kff_out
                    ? kff_writer.write(table, cfg)
                    : out ? write_counts<k, m>(table, cfg, chunk, *out, out_mutex)
                          : count_filtered<k, m>(table, cfg);
                total_written.fetch_add(wrt, std::memory_order_relaxed);

                // Collect per-partition overflow stats.
                const uint64_t ov_cnt = table.overflow_insert_count();
                if(ov_cnt > 0)
                {
                    ov_total.fetch_add(ov_cnt, std::memory_order_relaxed);
                    auto top = table.overflow_top_minimizers(20);
                    std::lock_guard<std::mutex> lg(ov_stats_mutex);
                    for(auto& [bin, cnt] : top)
                    {
                        auto it = std::find_if(ov_top_global.begin(), ov_top_global.end(),
                                               [bin](const auto& e){ return e.first == bin; });
                        if(it != ov_top_global.end()) it->second += cnt;
                        else ov_top_global.emplace_back(bin, cnt);
                    }
                }
            }
            kff_writer.flush();
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

    // Print overflow summary (always printed to stderr, regardless of -hp).
    if(ov_total.load() > 0)
    {
        const uint64_t tot_ov = ov_total.load();
        const uint64_t tot_ins = total_inserted.load();
        std::cerr << "[overflow] " << tot_ov << " k-mers went to overflow ("
                  << std::fixed << std::setprecision(2)
                  << (100.0 * tot_ov / tot_ins) << "% of total)\n";

        // Sort global top and print top-20.
        std::partial_sort(ov_top_global.begin(),
                          ov_top_global.begin() + std::min(std::size_t(20), ov_top_global.size()),
                          ov_top_global.end(),
                          [](const auto& a, const auto& b){ return a.second > b.second; });
        if(ov_top_global.size() > 20) ov_top_global.resize(20);

        std::cerr << "[overflow] top routing-hash bins (top "
                  << table_t::OV_HIST_BITS_PUBLIC << " bits):\n";
        std::cerr << "  rank     bin_id (hex)     overflow_kmers\n";
        for(std::size_t i = 0; i < ov_top_global.size(); ++i)
            std::cerr << "  " << std::setw(4) << (i+1)
                      << "     0x" << std::hex << std::setw(4) << std::setfill('0') << ov_top_global[i].first
                      << std::dec << std::setfill(' ')
                      << "     " << ov_top_global[i].second << "\n";
    }

    if (cfg.debug_stats)
        emit_debug_stats(part_infos, n_parts, cfg);
    if (cfg.debug_stats && cfg.dedup_mode != SuperkmerDedupMode::Off)
        emit_superkmer_dedup_stats(dedup_stats);

    return { total_inserted.load(), total_written.load() };
}


// ─── In-memory counting harness ───────────────────────────────────────────────
//
// Same as count_and_write but reads from per-partition std::string buffers
// (populated by partition_kmers_mem) instead of mmap'd disk files.
// Each partition buffer is cleared after processing to release memory
// incrementally — peak RSS ≈ largest-single-partition buffer, not all at once.

template <uint16_t k, uint16_t m, bool canonical_ = true>
std::pair<uint64_t, uint64_t> count_and_write_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    std::ofstream*            out,       // non-null for TSV output
    KffOutput*                kff_out)   // non-null for KFF output
{
    using table_t = phase2_table_t<k, false, uint32_t, m, canonical_>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    const size_t dedup_budget = superkmer_dedup_worker_budget(cfg, n_threads);
    std::atomic<size_t> next_part{0};

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<uint64_t> ov_total{0};

    // Debug statistics (only used when cfg.debug_stats).
    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};
    SuperkmerDedupStats dedup_stats;

    auto worker = [&](size_t /*tid*/) {
        typename table_t::Token token;
        std::string chunk;
        KffBatchWriter<k> kff_writer(kff_out);

        while (true) {
            const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
            if (p >= n_parts) break;
            const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
            size_t init_sz;
            if (cal > 0) {
                init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                     size_t(1u << 16), size_t(1u << 22));
            } else {
                const size_t per_part = (total_kmers > 0)
                    ? static_cast<size_t>(total_kmers / n_parts) * 2
                    : size_t(1u << 27) / n_parts;
                init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
            }
            table_t table(init_sz, 1);

            PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;

            uint64_t ins;
            {
                MemoryReader<k, m> reader(part_bufs[p]);
                ins = cfg.dedup_mode == SuperkmerDedupMode::Off
                    ? count_partition<k, m, canonical_, MemoryReader<k, m>>(
                        reader, table, token, dbg)
                    : count_partition_aggregated<k, m, canonical_>(
                        reader, cfg.dedup_mode, dedup_budget,
                        table, token, dbg,
                        cfg.debug_stats ? &dedup_stats : nullptr);
            }
            // Release the buffer immediately after counting to cap peak RSS.
            { std::string tmp; part_bufs[p].swap(tmp); }

            total_inserted.fetch_add(ins, std::memory_order_relaxed);
            ov_total.fetch_add(table.overflow_insert_count(), std::memory_order_relaxed);
            if (cal == 0) {
                const uint64_t unique = static_cast<uint64_t>(table.size());
                if (unique > 0) {
                    uint64_t expected = 0;
                    calibrated_unique.compare_exchange_strong(
                        expected, unique,
                        std::memory_order_relaxed, std::memory_order_relaxed);
                }
            }

            if (dbg) {
                dbg->init_sz     = init_sz;
                dbg->n_inserted  = ins;
                dbg->n_unique    = static_cast<uint64_t>(table.size());
                dbg->n_overflow  = table.overflow_insert_count();
                dbg->table_cap   = table.capacity();
                dbg->n_buckets   = table.bucket_count();
                dbg->resize_log  = table.resize_log();
            }

            const uint64_t wrt = kff_out
                ? kff_writer.write(table, cfg)
                : out ? write_counts<k, m>(table, cfg, chunk, *out, out_mutex)
                      : count_filtered<k, m>(table, cfg);
            total_written.fetch_add(wrt, std::memory_order_relaxed);
        }
        kff_writer.flush();
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    if (ov_total.load() > 0)
        std::cerr << "[overflow] " << ov_total.load() << " k-mers went to overflow\n";

    if (cfg.debug_stats)
        emit_debug_stats(part_infos, n_parts, cfg);
    if (cfg.debug_stats && cfg.dedup_mode != SuperkmerDedupMode::Off)
        emit_superkmer_dedup_stats(dedup_stats);

    return { total_inserted.load(), total_written.load() };
}
