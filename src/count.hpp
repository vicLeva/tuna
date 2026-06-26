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
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <exception>
#include <stdexcept>
#include <charconv>
#include <string_view>
#include <filesystem>
#include <chrono>

#include <ankerl/unordered_dense.h>
#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"


// Returns the smallest power of two >= v (returns 1 for v == 0).
inline size_t next_pow2(size_t v) noexcept {
    size_t p = 1;
    while (p < v) p <<= 1;
    return p;
}

inline uint64_t elapsed_ns(std::chrono::steady_clock::time_point t0) noexcept
{
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0).count());
}

inline void atomic_fetch_max_relaxed(std::atomic<uint64_t>& dst, uint64_t value) noexcept
{
    uint64_t cur = dst.load(std::memory_order_relaxed);
    while (cur < value &&
           !dst.compare_exchange_weak(
               cur, value, std::memory_order_relaxed, std::memory_order_relaxed)) {}
}

inline double ns_to_s(uint64_t ns) noexcept
{
    return static_cast<double>(ns) / 1'000'000'000.0;
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

struct DiskDedupStats {
    std::atomic<uint64_t> aggregate_attempts{0};
    std::atomic<uint64_t> aggregate_successes{0};
    std::atomic<uint64_t> aggregate_fallbacks{0};
    std::atomic<uint64_t> aggregate_records{0};
    std::atomic<uint64_t> aggregate_unique_records{0};
    std::atomic<uint64_t> aggregate_read_ns{0};
    std::atomic<uint64_t> aggregate_replay_ns{0};
    std::atomic<uint64_t> shard_rounds{0};
    std::atomic<uint64_t> shard_files{0};
    std::atomic<uint64_t> shard_input_bytes{0};
    std::atomic<uint64_t> shard_output_bytes{0};
};

struct Phase2TimingStats {
    std::atomic<uint64_t> partitions{0};
    std::atomic<uint64_t> partition_bytes{0};
    std::atomic<uint64_t> table_init_ns{0};
    std::atomic<uint64_t> count_ns{0};
    std::atomic<uint64_t> output_ns{0};
    std::atomic<uint64_t> max_part_ns{0};
    std::atomic<uint64_t> max_count_ns{0};
    std::atomic<uint64_t> max_output_ns{0};
};


// ─── Counting brick ───────────────────────────────────────────────────────────
//
// Drains a superkmer reader into the caller-provided hash table (mt_=false).
// The stored min_pos byte lets Phase 2 init the bucket hash in O(m) rather than
// O(k).  After each superkmer window is initialized, we prefetch the destination
// bucket from the already-computed minimizer hash.
//
// Returns the total number of k-mer insertions (with multiplicity).

template <uint16_t k, uint16_t m, typename Reader = SuperkmerReader<k, m>>
uint64_t count_partition(
    Reader&                                                               reader,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token,
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

    kache_hash::Kmer_Window<k, m> win;
    if (cur_len >= k) {
        if (cur_min_pos != NO_MIN) {
            const uint64_t mh = win.init_packed_with_min(cur_packed, cur_min_pos);
            if (dbg) min_kmer_count[mh] += static_cast<uint32_t>(cur_len - k + 1);
        } else {
            win.init_packed(cur_packed);
        }
        table.prefetch(win);  // cold-start: no previous superkmer to hide behind
    }

    // ── Main loop ─────────────────────────────────────────────────────────────
    while (reader.next()) {
        const uint8_t* nxt_packed  = reader.packed_data();  // packed bases of next superkmer (prefetch target)
        const size_t   nxt_len     = reader.size();
        const hdr_t    nxt_min_pos = reader.min_pos();

        // Process CURRENT superkmer.
        if (cur_len >= k) {
            table.upsert(win, inc, uint32_t(1), token);
            ++inserted;

            // Unpack subsequent bases directly as DNA::Base (kache encoding).
            const uint8_t* byte_ptr = cur_packed + (k >> 2);
            int shift = static_cast<int>(6u - 2u * (k & 3u));
            for (size_t i = k; i < cur_len; ++i) {
                const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
                shift -= 2;
                if (shift < 0) { shift = 6; ++byte_ptr; }
                win.advance(b);
                table.upsert(win, inc, uint32_t(1), token);
                ++inserted;
            }
        }

        // Advance to next: reinitialise window from the already-prefetched data.
        cur_packed  = nxt_packed;
        cur_len     = nxt_len;
        cur_min_pos = nxt_min_pos;
        if (cur_len >= k) {
            if (cur_min_pos != NO_MIN) {
                const uint64_t mh = win.init_packed_with_min(cur_packed, cur_min_pos);
                if (dbg) min_kmer_count[mh] += static_cast<uint32_t>(cur_len - k + 1);
                table.prefetch(win);
            } else {
                win.init_packed(cur_packed);
                table.prefetch(win);  // NO_MIN fallback: prefetch after init
            }
        }
    }

    // ── Last superkmer ────────────────────────────────────────────────────────
    if (cur_len >= k) {
        table.upsert(win, inc, uint32_t(1), token);
        ++inserted;
        const uint8_t* byte_ptr = cur_packed + (k >> 2);
        int shift = static_cast<int>(6u - 2u * (k & 3u));
        for (size_t i = k; i < cur_len; ++i) {
            const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
            shift -= 2;
            if (shift < 0) { shift = 6; ++byte_ptr; }
            win.advance(b);
            table.upsert(win, inc, uint32_t(1), token);
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


template <uint16_t k, uint16_t m>
uint64_t count_superkmer_record(
    const uint8_t*                                                        packed,
    size_t                                                                len,
    sk_hdr_t<k, m>                                                        min_pos,
    uint32_t                                                              multiplicity,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token)
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    if (len < k || multiplicity == 0) return 0;

    auto inc = [multiplicity](uint32_t v) { return v + multiplicity; };
    uint64_t inserted = 0;

    kache_hash::Kmer_Window<k, m> win;
    if (min_pos != NO_MIN) {
        win.init_packed_with_min(packed, min_pos);
        table.prefetch(win);
    } else {
        win.init_packed(packed);
        table.prefetch(win);
    }

    table.upsert(win, inc, multiplicity, token);
    inserted += multiplicity;

    const uint8_t* byte_ptr = packed + (k >> 2);
    int shift = static_cast<int>(6u - 2u * (k & 3u));
    for (size_t i = k; i < len; ++i) {
        const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
        shift -= 2;
        if (shift < 0) { shift = 6; ++byte_ptr; }
        win.advance(b);
        table.upsert(win, inc, multiplicity, token);
        inserted += multiplicity;
    }

    return inserted;
}


template <uint16_t k, uint16_t m>
struct PackedSuperkmerKey {
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);
    static constexpr size_t MAX_LEN = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr size_t MAX_RECORD_BYTES = HDR_BYTES + (MAX_LEN + 3u) / 4u;
    static constexpr size_t WORDS = (MAX_RECORD_BYTES + sizeof(uint64_t) - 1u) / sizeof(uint64_t);
    static constexpr size_t STORAGE_BYTES = WORDS * sizeof(uint64_t);

    std::array<uint64_t, WORDS> words{};

    static PackedSuperkmerKey from_record(const char* record, size_t record_size)
    {
        PackedSuperkmerKey key;
        if (record_size > MAX_RECORD_BYTES)
            throw std::runtime_error("tuna: corrupt superkmer record exceeds maximum encoded length");
        std::memcpy(key.words.data(), record, record_size);
        return key;
    }

    const char* bytes() const noexcept
    {
        return reinterpret_cast<const char*>(words.data());
    }

    friend bool operator==(const PackedSuperkmerKey& a,
                           const PackedSuperkmerKey& b) noexcept
    {
        return a.words == b.words;
    }
};


template <uint16_t k, uint16_t m>
struct PackedSuperkmerKeyHash {
    size_t operator()(const PackedSuperkmerKey<k, m>& key) const noexcept
    {
        return static_cast<size_t>(
            XXH3_64bits(key.words.data(), PackedSuperkmerKey<k, m>::STORAGE_BYTES));
    }
};


template <uint16_t k, uint16_t m>
uint64_t count_packed_superkmer_key(
    const PackedSuperkmerKey<k, m>&                                      key,
    uint32_t                                                             multiplicity,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&        table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token)
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    hdr_t len, mp;
    std::memcpy(&len, key.bytes(),                 sizeof(hdr_t));
    std::memcpy(&mp,  key.bytes() + sizeof(hdr_t), sizeof(hdr_t));
    const auto* packed = reinterpret_cast<const uint8_t*>(key.bytes() + HDR_BYTES);
    return count_superkmer_record<k, m>(
        packed, static_cast<size_t>(len), mp, multiplicity, table, token);
}


inline uint64_t superkmer_file_size(const std::string& path)
{
    std::error_code ec;
    const auto sz = std::filesystem::file_size(path, ec);
    return ec ? 0 : static_cast<uint64_t>(sz);
}


inline std::string dedup_shard_path(
    const std::string& work_dir,
    size_t             partition,
    size_t             level,
    size_t             shard)
{
    return work_dir + "dedup_p" + std::to_string(partition)
        + "_l" + std::to_string(level)
        + "_s" + std::to_string(shard) + ".superkmers";
}


inline size_t disk_dedup_worker_budget(const Config& cfg, size_t n_threads)
{
    static constexpr uint64_t MiB = uint64_t(1) << 20;
    if (cfg.ram_budget_bytes == 0)
        return static_cast<size_t>(256 * MiB);

    const uint64_t threads = std::max<uint64_t>(1, n_threads);
    // This budget is only for the exact superkmer dedup map.  The phase-2
    // k-mer table and other working memory are separate, so keep substantial
    // headroom while still using enough of -ram to avoid excess sharding.
    const uint64_t share = (cfg.ram_budget_bytes * 3u / 8u) / threads;
    return static_cast<size_t>(std::min(std::max<uint64_t>(share, 16 * MiB), 512 * MiB));
}


template <uint16_t k, uint16_t m>
bool try_count_disk_aggregated(
    const std::string&                                                    path,
    size_t                                                               memory_budget,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token,
    uint64_t&                                                            inserted,
    DiskDedupStats*                                                       stats = nullptr)
{
    using key_t = PackedSuperkmerKey<k, m>;
    using map_t = ankerl::unordered_dense::map<key_t, uint32_t, PackedSuperkmerKeyHash<k, m>>;
    static constexpr size_t ENTRY_OVERHEAD = 32;
    static constexpr size_t ENTRY_ESTIMATE = sizeof(key_t) + sizeof(uint32_t) + ENTRY_OVERHEAD;

    inserted = 0;
    if (superkmer_file_size(path) == 0) return true;

    SuperkmerReader<k, m> reader(path);
    if (!reader.ok())
        throw std::runtime_error("tuna: cannot open partition file for reading: " + path);

    map_t multiplicities;
    if (memory_budget > ENTRY_ESTIMATE) {
        const size_t reserve_by_budget = memory_budget / ENTRY_ESTIMATE;
        const size_t reserve_by_file = static_cast<size_t>(
            std::max<uint64_t>(1024, superkmer_file_size(path) / key_t::MAX_RECORD_BYTES));
        multiplicities.reserve(std::min(reserve_by_budget, reserve_by_file));
    }

    size_t total_key_bytes = 0;
    size_t records_seen = 0;
    const auto read_t0 = std::chrono::steady_clock::now();
    while (reader.next()) {
        const auto key = key_t::from_record(reader.record_data(), reader.record_size());
        auto [it, is_new] = multiplicities.try_emplace(key, uint32_t(1));
        if (is_new) {
            total_key_bytes += key_t::STORAGE_BYTES;
        } else {
            ++it->second;
        }

        if ((++records_seen & 4095u) == 0) {
            const size_t estimated = multiplicities.size() * ENTRY_OVERHEAD
                                   + total_key_bytes
                                   + multiplicities.size() * sizeof(uint32_t);
            if (estimated > memory_budget) return false;
        }
    }
    const uint64_t read_ns = elapsed_ns(read_t0);

    const auto replay_t0 = std::chrono::steady_clock::now();
    for (const auto& [key, multiplicity] : multiplicities)
        inserted += count_packed_superkmer_key<k, m>(key, multiplicity, table, token);
    const uint64_t replay_ns = elapsed_ns(replay_t0);

    if (stats) {
        stats->aggregate_records.fetch_add(records_seen, std::memory_order_relaxed);
        stats->aggregate_unique_records.fetch_add(multiplicities.size(), std::memory_order_relaxed);
        stats->aggregate_read_ns.fetch_add(read_ns, std::memory_order_relaxed);
        stats->aggregate_replay_ns.fetch_add(replay_ns, std::memory_order_relaxed);
    }

    return true;
}


template <uint16_t k, uint16_t m>
uint64_t count_disk_dedup_exact(
    const std::string&                                                    path,
    const Config&                                                         cfg,
    size_t                                                               partition,
    size_t                                                               level,
    size_t                                                               memory_budget,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token,
    DiskDedupStats*                                                       stats = nullptr)
{
    static constexpr uint64_t MEMORY_FACTOR = 4;
    static constexpr size_t MAX_SHARDS = 64;

    const uint64_t file_bytes = superkmer_file_size(path);
    if (file_bytes == 0) return 0;

    uint64_t inserted = 0;
    if (stats) stats->aggregate_attempts.fetch_add(1, std::memory_order_relaxed);
    if (try_count_disk_aggregated<k, m>(path, memory_budget, table, token, inserted, stats)) {
        if (stats) stats->aggregate_successes.fetch_add(1, std::memory_order_relaxed);
        return inserted;
    }
    if (stats) stats->aggregate_fallbacks.fetch_add(1, std::memory_order_relaxed);

    if (level >= 12)
        throw std::runtime_error("tuna: disk dedup shard exceeded memory budget after recursive sharding");

    const uint64_t target = std::max<uint64_t>(1, memory_budget);
    const uint64_t needed = std::max<uint64_t>(2, (file_bytes * MEMORY_FACTOR + target - 1) / target);
    const size_t n_shards = std::min(MAX_SHARDS, next_pow2(static_cast<size_t>(needed)));
    if (stats) {
        stats->shard_rounds.fetch_add(1, std::memory_order_relaxed);
        stats->shard_files.fetch_add(n_shards, std::memory_order_relaxed);
        stats->shard_input_bytes.fetch_add(file_bytes, std::memory_order_relaxed);
    }

    std::vector<std::string> shard_paths;
    shard_paths.reserve(n_shards);
    for (size_t s = 0; s < n_shards; ++s)
        shard_paths.push_back(dedup_shard_path(cfg.work_dir, partition, level, s));

    auto cleanup = [&]() {
        for (const auto& shard : shard_paths) {
            std::error_code ec;
            std::filesystem::remove(shard, ec);
        }
    };

    try {
        std::vector<std::ofstream> outs;
        outs.reserve(n_shards);
        for (const auto& shard : shard_paths) {
            outs.emplace_back(shard, std::ios::binary | std::ios::trunc);
            if (!outs.back())
                throw std::runtime_error("tuna: cannot open dedup shard for writing: " + shard);
        }

        SuperkmerReader<k, m> reader(path);
        if (!reader.ok())
            throw std::runtime_error("tuna: cannot open partition file for reading: " + path);

#ifdef TUNA_LZ4_BUCKETS
        static constexpr size_t SHARD_BUFFER_BYTES = 1u << 20;
        std::vector<std::vector<char>> shard_buffers(n_shards);
        std::vector<uint8_t> shard_initialized(n_shards, 0);
        auto flush_shard = [&](size_t shard) {
            auto& buf = shard_buffers[shard];
            if (buf.empty()) return;
            if (!shard_initialized[shard]) {
                tuna_superkmer_detail::write_lz4_bucket_magic(outs[shard]);
                shard_initialized[shard] = 1;
            }
            tuna_superkmer_detail::write_superkmer_bucket_block(
                outs[shard], buf.data(), buf.size(), true);
            buf.clear();
        };
#endif

        const uint64_t seed = 0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(level);
        while (reader.next()) {
            const uint64_t h = XXH3_64bits_withSeed(
                reader.record_data(), reader.record_size(), seed);
            const size_t shard = static_cast<size_t>(h) & (n_shards - 1u);
#ifdef TUNA_LZ4_BUCKETS
            if (cfg.lz4_shards) {
                auto& buf = shard_buffers[shard];
                buf.insert(buf.end(), reader.record_data(), reader.record_data() + reader.record_size());
                if (buf.size() >= SHARD_BUFFER_BYTES)
                    flush_shard(shard);
            } else
#endif
            outs[shard].write(reader.record_data(),
                              static_cast<std::streamsize>(reader.record_size()));
        }
#ifdef TUNA_LZ4_BUCKETS
        if (cfg.lz4_shards) {
            for (size_t s = 0; s < n_shards; ++s)
                flush_shard(s);
        }
#endif
        for (auto& out : outs) {
            out.close();
            if (!out)
                throw std::runtime_error("tuna: failed while writing dedup shard");
        }
        if (stats) {
            uint64_t written = 0;
            for (const auto& shard : shard_paths) {
                std::error_code ec;
                const auto sz = std::filesystem::file_size(shard, ec);
                if (!ec) written += static_cast<uint64_t>(sz);
            }
            stats->shard_output_bytes.fetch_add(written, std::memory_order_relaxed);
        }

        uint64_t total_inserted = 0;
        for (size_t s = 0; s < n_shards; ++s) {
            total_inserted += count_disk_dedup_exact<k, m>(
                shard_paths[s], cfg, partition, level + 1, memory_budget, table, token, stats);
            std::error_code ec;
            std::filesystem::remove(shard_paths[s], ec);
        }
        return total_inserted;
    } catch (...) {
        cleanup();
        throw;
    }
}


template <uint16_t k, uint16_t m>
uint64_t count_partition_mem_aggregated(
    const std::string&                                                    data,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>::Token& token)
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    struct RecordHash {
        size_t operator()(std::string_view s) const noexcept {
            return static_cast<size_t>(XXH3_64bits(s.data(), s.size()));
        }
    };

    ankerl::unordered_dense::map<std::string_view, uint32_t, RecordHash> multiplicities;
    multiplicities.reserve(std::max<size_t>(1024, data.size() / 8));

    const char* cur = data.data();
    const char* end = data.data() + data.size();
    while (cur + static_cast<ptrdiff_t>(HDR_BYTES) <= end) {
        const char* rec = cur;
        hdr_t len, mp;
        std::memcpy(&len, cur,                 sizeof(hdr_t));
        std::memcpy(&mp,  cur + sizeof(hdr_t), sizeof(hdr_t));
        (void)mp;
        if (len == 0) break;
        cur += HDR_BYTES;
        const size_t packed_bytes = (static_cast<size_t>(len) + 3u) / 4u;
        if (cur + static_cast<ptrdiff_t>(packed_bytes) > end) break;
        cur += packed_bytes;
        ++multiplicities[std::string_view(rec, static_cast<size_t>(cur - rec))];
    }

    uint64_t inserted = 0;
    for (const auto& [record, multiplicity] : multiplicities) {
        hdr_t len, mp;
        std::memcpy(&len, record.data(),                 sizeof(hdr_t));
        std::memcpy(&mp,  record.data() + sizeof(hdr_t), sizeof(hdr_t));
        const auto* packed = reinterpret_cast<const uint8_t*>(record.data() + HDR_BYTES);
        inserted += count_superkmer_record<k, m>(
            packed, static_cast<size_t>(len), mp, multiplicity, table, token);
    }

    return inserted;
}


// ─── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m, bool mt_ = false>
uint64_t write_counts(
    kache_hash::Streaming_Kmer_Hash_Table<k, mt_, uint32_t, m>& table,
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

template <uint16_t k, uint16_t m, bool mt_ = false>
uint64_t count_output_kmers(
    kache_hash::Streaming_Kmer_Hash_Table<k, mt_, uint32_t, m>& table,
    const Config& cfg)
{
    if (cfg.ci <= 1 && cfg.cx == std::numeric_limits<uint64_t>::max())
        return static_cast<uint64_t>(table.size());

    uint64_t written = 0;
    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        if (cnt >= cfg.ci && cnt <= cfg.cx) ++written;
    });
    return written;
}


// ─── KFF output brick ─────────────────────────────────────────────────────────
//
// Encodes each k-mer from the table as 2-bit MSB-first bytes and flushes to
// a KffOutput in batches of ~1 MB.  Thread-safe via KffOutput::write_batch.

template <uint16_t k, uint16_t m, bool mt_ = false>
uint64_t write_counts_kff(
    kache_hash::Streaming_Kmer_Hash_Table<k, mt_, uint32_t, m>& table,
    const Config& cfg,
    KffOutput&    kff_out)
{
    constexpr size_t KMER_BYTES  = (k + 3) / 4;                            // bytes per k-mer (2-bit packed)
    constexpr size_t BATCH_KMERS = (size_t(1) << 20) / (KMER_BYTES + 4);  // k-mers per KFF write batch

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

    table.for_each([&](const auto& entry) {
        const uint64_t cnt = entry.second;
        if (cnt < cfg.ci || cnt > cfg.cx) return;

        // Pack label into 2-bit bytes (MSB = first base).
        seq_buf.resize(seq_buf.size() + KMER_BYTES, 0);
        uint8_t* dst = seq_buf.data() + seq_buf.size() - KMER_BYTES;
        entry.first.write_packed_2bit_msb(dst);

        cnt_buf.push_back(static_cast<uint32_t>(cnt));
        ++written;

        if (cnt_buf.size() >= BATCH_KMERS) flush();
    });
    flush();

    return written;
}


// ─── Callback output brick ────────────────────────────────────────────────────
//
// Drains the table, applies ci/cx filters, and calls cb(kmer, count) for each
// passing k-mer.  cb is called from the calling thread only (one table per call),
// so no internal mutex is used.  If multiple worker threads drain different tables
// concurrently, cb may be invoked from several threads simultaneously — the caller
// is responsible for any needed synchronisation.

template <uint16_t k, uint16_t m, bool mt_ = false, typename Callback>
uint64_t write_counts_callback(
    kache_hash::Streaming_Kmer_Hash_Table<k, mt_, uint32_t, m>& table,
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

template <uint16_t k, uint16_t m, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    Callback&&                cb)
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>;  // hash table type (local alias)

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
                                         size_t(1u << 15), size_t(1u << 22));
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
                    ins = count_partition<k, m, MemoryReader<k, m>>(reader, table, token);
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


template <uint16_t k, uint16_t m, typename Callback>
std::pair<uint64_t, uint64_t> count_and_callback(
    const Config& cfg,
    uint64_t      total_kmers,
    Callback&&    cb)
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    const size_t dedup_budget = disk_dedup_worker_budget(cfg, n_threads);
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
                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                         size_t(1u << 15), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }
                table_t table(init_sz, 1);

                const uint64_t ins = count_disk_dedup_exact<k, m>(
                    path, cfg, p, 0, dedup_budget, table, token);
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

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    uint64_t       total_kmers,
    std::ofstream* out,       // non-null for TSV output
    KffOutput*     kff_out)   // non-null for KFF output
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    const size_t dedup_budget = disk_dedup_worker_budget(cfg, n_threads);
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
    const bool collect_phase2_stats = cfg.debug_stats;
    Phase2TimingStats phase2_stats;
    DiskDedupStats dedup_stats;

    auto worker = [&](size_t /*tid*/) {
        try {
            typename table_t::Token token;
            std::string chunk;

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
                if (p >= n_parts) break;
                const std::string path = partition_path(cfg.work_dir, p);
                const uint64_t part_bytes = collect_phase2_stats ? superkmer_file_size(path) : 0;
                const auto part_t0 = collect_phase2_stats
                    ? std::chrono::steady_clock::now()
                    : std::chrono::steady_clock::time_point{};
                const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
                size_t init_sz;
                if (cal > 0) {
                    init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                         size_t(1u << 15), size_t(1u << 22));
                } else {
                    const size_t per_part = (total_kmers > 0)
                        ? static_cast<size_t>(total_kmers / n_parts) * 2
                        : size_t(1u << 27) / n_parts;
                    init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
                }
                const auto table_t0 = collect_phase2_stats
                    ? std::chrono::steady_clock::now()
                    : std::chrono::steady_clock::time_point{};
                table_t table(init_sz, 1);
                const uint64_t table_ns = collect_phase2_stats ? elapsed_ns(table_t0) : 0;

                PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;
                uint64_t ins;
                const auto count_t0 = collect_phase2_stats
                    ? std::chrono::steady_clock::now()
                    : std::chrono::steady_clock::time_point{};
                if (!cfg.debug_stats) {
                    ins = count_disk_dedup_exact<k, m>(
                        path, cfg, p, 0, dedup_budget, table, token,
                        collect_phase2_stats ? &dedup_stats : nullptr);
                } else {
                    SuperkmerReader<k, m> reader(path);
                    if (!reader.ok())
                        throw std::runtime_error("tuna: cannot open partition file for reading: " + path);
                    ins = count_partition<k, m>(reader, table, token, dbg);
                }
                const uint64_t count_ns = collect_phase2_stats ? elapsed_ns(count_t0) : 0;
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

                const auto output_t0 = collect_phase2_stats
                    ? std::chrono::steady_clock::now()
                    : std::chrono::steady_clock::time_point{};
                const uint64_t wrt = cfg.count_only
                    ? count_output_kmers<k, m>(table, cfg)
                    : (kff_out
                        ? write_counts_kff<k, m>(table, cfg, *kff_out)
                        : write_counts<k, m>(table, cfg, chunk, *out, out_mutex));
                const uint64_t output_ns = collect_phase2_stats ? elapsed_ns(output_t0) : 0;
                total_written.fetch_add(wrt, std::memory_order_relaxed);
                if (collect_phase2_stats) {
                    const uint64_t part_ns = elapsed_ns(part_t0);
                    phase2_stats.partitions.fetch_add(1, std::memory_order_relaxed);
                    phase2_stats.partition_bytes.fetch_add(part_bytes, std::memory_order_relaxed);
                    phase2_stats.table_init_ns.fetch_add(table_ns, std::memory_order_relaxed);
                    phase2_stats.count_ns.fetch_add(count_ns, std::memory_order_relaxed);
                    phase2_stats.output_ns.fetch_add(output_ns, std::memory_order_relaxed);
                    atomic_fetch_max_relaxed(phase2_stats.max_part_ns, part_ns);
                    atomic_fetch_max_relaxed(phase2_stats.max_count_ns, count_ns);
                    atomic_fetch_max_relaxed(phase2_stats.max_output_ns, output_ns);
                }

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

    if (collect_phase2_stats) {
        const uint64_t p_done = phase2_stats.partitions.load(std::memory_order_relaxed);
        const double mb = static_cast<double>(
            phase2_stats.partition_bytes.load(std::memory_order_relaxed)) / (1024.0 * 1024.0);
        std::cerr << "[phase2] partitions=" << p_done
                  << " input_mib=" << std::fixed << std::setprecision(1) << mb
                  << " table_init=" << std::setprecision(3)
                  << ns_to_s(phase2_stats.table_init_ns.load(std::memory_order_relaxed)) << "s"
                  << " count=" << ns_to_s(phase2_stats.count_ns.load(std::memory_order_relaxed)) << "s"
                  << " output_scan=" << ns_to_s(phase2_stats.output_ns.load(std::memory_order_relaxed)) << "s"
                  << " max_part=" << ns_to_s(phase2_stats.max_part_ns.load(std::memory_order_relaxed)) << "s"
                  << " max_count=" << ns_to_s(phase2_stats.max_count_ns.load(std::memory_order_relaxed)) << "s"
                  << " max_output=" << ns_to_s(phase2_stats.max_output_ns.load(std::memory_order_relaxed)) << "s"
                  << "\n";

        std::cerr << "[phase2-dedup] budget_mib="
                  << (static_cast<double>(dedup_budget) / (1024.0 * 1024.0))
                  << " attempts=" << dedup_stats.aggregate_attempts.load(std::memory_order_relaxed)
                  << " direct=" << dedup_stats.aggregate_successes.load(std::memory_order_relaxed)
                  << " sharded=" << dedup_stats.aggregate_fallbacks.load(std::memory_order_relaxed)
                  << " shard_rounds=" << dedup_stats.shard_rounds.load(std::memory_order_relaxed)
                  << " shard_files=" << dedup_stats.shard_files.load(std::memory_order_relaxed)
                  << " shard_input_mib=" << std::setprecision(1)
                  << (static_cast<double>(dedup_stats.shard_input_bytes.load(std::memory_order_relaxed)) /
                      (1024.0 * 1024.0))
                  << " shard_output_mib="
                  << (static_cast<double>(dedup_stats.shard_output_bytes.load(std::memory_order_relaxed)) /
                      (1024.0 * 1024.0))
                  << "\n";
        const uint64_t records = dedup_stats.aggregate_records.load(std::memory_order_relaxed);
        if (records > 0) {
            const uint64_t unique_records = dedup_stats.aggregate_unique_records.load(std::memory_order_relaxed);
            std::cerr << "[phase2-dedup-detail] records=" << records
                      << " unique_records=" << unique_records
                      << " multiplicity=" << std::setprecision(3)
                      << (static_cast<double>(records) / static_cast<double>(std::max<uint64_t>(1, unique_records)))
                      << " aggregate_read=" << ns_to_s(dedup_stats.aggregate_read_ns.load(std::memory_order_relaxed)) << "s"
                      << " aggregate_replay=" << ns_to_s(dedup_stats.aggregate_replay_ns.load(std::memory_order_relaxed)) << "s"
                      << "\n";
        }
    }

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

        std::cerr << "[overflow] top minimizer bins (top " << table_t::OV_HIST_BITS_PUBLIC
                  << " bits of ntHash canonical):\n";
        std::cerr << "  rank     bin_id (hex)     overflow_kmers\n";
        for(std::size_t i = 0; i < ov_top_global.size(); ++i)
            std::cerr << "  " << std::setw(4) << (i+1)
                      << "     0x" << std::hex << std::setw(4) << std::setfill('0') << ov_top_global[i].first
                      << std::dec << std::setfill(' ')
                      << "     " << ov_top_global[i].second << "\n";
    }

    if (cfg.debug_stats)
        emit_debug_stats(part_infos, n_parts, cfg);

    return { total_inserted.load(), total_written.load() };
}


// ─── In-memory counting harness ───────────────────────────────────────────────
//
// Same as count_and_write but reads from per-partition std::string buffers
// (populated by partition_kmers_mem) instead of mmap'd disk files.
// Each partition buffer is cleared after processing to release memory
// incrementally — peak RSS ≈ largest-single-partition buffer, not all at once.

template <uint16_t k, uint16_t m>
std::pair<uint64_t, uint64_t> count_and_write_mem(
    const Config&             cfg,
    uint64_t                  total_kmers,
    std::vector<std::string>& part_bufs,
    std::ofstream*            out,       // non-null for TSV output
    KffOutput*                kff_out)   // non-null for KFF output
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, m>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);
    std::atomic<size_t> next_part{0};

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};
    std::atomic<uint64_t> ov_total{0};

    // Debug statistics (only used when cfg.debug_stats).
    std::vector<PartitionDebugInfo> part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);
    std::atomic<uint64_t> calibrated_unique{0};

    auto worker = [&](size_t /*tid*/) {
        typename table_t::Token token;
        std::string chunk;

        while (true) {
            const size_t p = next_part.fetch_add(1, std::memory_order_relaxed);
            if (p >= n_parts) break;
            const uint64_t cal = calibrated_unique.load(std::memory_order_relaxed);
            size_t init_sz;
            if (cal > 0) {
                init_sz = std::clamp(next_pow2(static_cast<size_t>(cal / 0.75) + 1),
                                     size_t(1u << 15), size_t(1u << 22));
            } else {
                const size_t per_part = (total_kmers > 0)
                    ? static_cast<size_t>(total_kmers / n_parts) * 2
                    : size_t(1u << 27) / n_parts;
                init_sz = std::clamp(per_part, size_t(1u << 15), size_t(1u << 22));
            }
            table_t table(init_sz, 1);

            PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;

            uint64_t ins;
            if (!cfg.debug_stats) {
                ins = count_partition_mem_aggregated<k, m>(part_bufs[p], table, token);
            } else {
                MemoryReader<k, m> reader(part_bufs[p]);
                ins = count_partition<k, m, MemoryReader<k, m>>(reader, table, token, dbg);
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

            const uint64_t wrt = cfg.count_only
                ? count_output_kmers<k, m>(table, cfg)
                : (kff_out
                    ? write_counts_kff<k, m>(table, cfg, *kff_out)
                    : write_counts<k, m>(table, cfg, chunk, *out, out_mutex));
            total_written.fetch_add(wrt, std::memory_order_relaxed);
        }
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

    return { total_inserted.load(), total_written.load() };
}
