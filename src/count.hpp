#pragma once

// Phase 2+3 — partition files → counted k-mers → output.
//
// Three independently modifiable bricks:
//
//   count_partition<k,l>   — read one partition file, fill a hash table.
//   write_counts<k,l>      — drain one table to the output stream.
//   count_and_write<k,l>   — parallel harness (threading / scheduling).

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "nt_hash.hpp"

#include <fstream>
#include <string>
#include <mutex>
#include <thread>
#include <atomic>
#include <utility>
#include <algorithm>
#include <iomanip>

#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"


// ─── Counting brick ───────────────────────────────────────────────────────────
//
// Drain an already-open SuperkmerReader into the caller-provided hash table.
// The table must be private to the calling thread (mt_=false, no locking).
//
// The stored min_pos header byte lets Phase 2 compute ntHash(minimizer) in O(l)
// instead of running MinimizerWindow::reset() in O(k).  All k-mers in the
// superkmer share the same minimizer, so init_packed_with_hash sets the bucket
// hash once and advance() skips nt_min entirely.  A single prefetch() hides
// the LLC miss (~200 ns) behind the O(l) lmer hash computation.
//
// Returns the total number of k-mer insertions (with multiplicity).

template <uint16_t k, uint16_t l>
uint64_t count_partition(
    SuperkmerReader&                                                      reader,
    kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>&         table,
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>::Token& token)
{
    // Compute canonical ntHash of the l-mer at position min_pos in packed data.
    // Decodes l bases (kache→ASCII) then runs nt_hash::Roller.  O(l) work —
    // replaces O(k) MinimizerWindow::reset() in init_packed.
    auto lmer_nt_hash = [](const uint8_t* packed, uint8_t pos) -> uint64_t {
        static constexpr char KACHE2ASCII[4] = {'A', 'C', 'G', 'T'};
        char buf[l];
        for (uint16_t i = 0; i < l; ++i) {
            const uint8_t p = static_cast<uint8_t>(pos + i);
            buf[i] = KACHE2ASCII[(packed[p >> 2] >> (6u - 2u * (p & 3u))) & 3u];
        }
        nt_hash::Roller<l> roller;
        roller.init(buf);
        return roller.canonical();
    };

    auto inc = [](uint32_t v) { return v + 1; };
    uint64_t inserted = 0;

    while (reader.next()) {
        const uint8_t* packed   = reader.packed_data();
        const size_t   len      = reader.size();
        const uint8_t  min_pos  = reader.min_pos();
        if (len < k) continue;

        // min_pos == 0xFF is a sentinel (KMC mode: multiple ntHash minimizers
        // possible within one superkmer) — fall back to MinimizerWindow::reset().
        // Otherwise use the stored position to compute ntHash in O(l) vs O(k).
        kache_hash::Kmer_Window<k, l> win;
        if (min_pos != 0xFF)
            win.init_packed_with_hash(packed, lmer_nt_hash(packed, min_pos));
        else
            win.init_packed(packed);

        // Prefetch the primary bucket before the first upsert.
        table.prefetch(win);

        table.upsert(win, inc, uint32_t(1), token);
        ++inserted;

        // Unpack subsequent bases directly as DNA::Base (kache encoding), no
        // ASCII round-trip.  Maintain byte pointer + shift counter to avoid
        // division/modulo in the hot loop.
        const uint8_t* byte_ptr = packed + (k >> 2);
        int shift = static_cast<int>(6u - 2u * (k & 3u));
        for (size_t i = k; i < len; ++i) {
            const auto b = static_cast<kache_hash::DNA::Base>((*byte_ptr >> shift) & 3u);
            shift -= 2;
            if (shift < 0) { shift = 6; ++byte_ptr; }
            win.advance(b);
            table.upsert(win, inc, uint32_t(1), token);
            ++inserted;
        }
    }
    return inserted;
}


// ─── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t l, bool mt_ = false>
uint64_t write_counts(
    kache_hash::Streaming_Kmer_Hash_Table<k, mt_, uint32_t, l>& table,
    const Config&   cfg,
    std::string&    chunk,
    std::ofstream&  out,
    std::mutex&     out_mutex)
{
    constexpr size_t WRITE_BATCH = 1u << 20; // 1 MB

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
        chunk += std::to_string(cnt);
        chunk += '\n';
        ++written;
        if (chunk.size() >= WRITE_BATCH) flush_chunk();
    });
    if (!chunk.empty()) flush_chunk();
    return written;
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// init_sz is fixed at 1<<22 (4M k-mer slots).  Rationale:
//   - For high-coverage data (e.g. 200 E. coli × 32 partitions → ~1.5M unique
//     k-mers/partition): 1.5M << 3.2M resize threshold → no resize, table is
//     8× smaller than the old file_size/4 heuristic → better L3 utilisation.
//   - For low-coverage data (e.g. human 1 file × 128 partitions → ~23M unique
//     k-mers/partition): table resizes 2-3× but amortised cost is small
//     compared to the counting time.

template <uint16_t k, uint16_t l>
std::pair<uint64_t, uint64_t> count_and_write(
    const Config&  cfg,
    std::ofstream& out)
{
    using table_t = kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>;

    const size_t n_parts   = cfg.num_partitions;
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_parts);

    std::mutex            out_mutex;
    std::atomic<uint64_t> total_inserted{0}, total_written{0};

    // Overflow statistics accumulated across all partitions.
    std::mutex                                     ov_stats_mutex;
    std::vector<std::pair<uint32_t, uint64_t>>     ov_top_global;   // merged top minimizers
    std::atomic<uint64_t>                          ov_total{0};

    auto worker = [&](size_t tid) {
        typename table_t::Token token;
        std::string chunk;

        for (size_t p = tid; p < n_parts; p += n_threads) {
            const std::string part_path = cfg.work_dir + cfg.partition_prefix()
                                          + "_" + std::to_string(p) + ".superkmers";
            SuperkmerReader reader(part_path);
            table_t table(1u << 22, 1); // 4M k-mers, 1 resize worker

            const uint64_t ins = count_partition<k, l>(reader, table, token);
            total_inserted.fetch_add(ins, std::memory_order_relaxed);

            const uint64_t wrt = write_counts<k, l>(table, cfg, chunk, out, out_mutex);
            total_written.fetch_add(wrt, std::memory_order_relaxed);

            // Collect per-partition overflow stats.
            const uint64_t ov_cnt = table.overflow_insert_count();
            if(ov_cnt > 0)
            {
                ov_total.fetch_add(ov_cnt, std::memory_order_relaxed);
                auto top = table.overflow_top_minimizers(20);
                std::lock_guard<std::mutex> lg(ov_stats_mutex);
                // Merge into global top list.
                for(auto& [bin, cnt] : top)
                {
                    auto it = std::find_if(ov_top_global.begin(), ov_top_global.end(),
                                           [bin](const auto& e){ return e.first == bin; });
                    if(it != ov_top_global.end()) it->second += cnt;
                    else ov_top_global.emplace_back(bin, cnt);
                }
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

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

        std::cerr << "[overflow] top minimizer bins (top 16 bits of ntHash canonical):\n";
        std::cerr << "  rank     bin_id (hex)     overflow_kmers\n";
        for(std::size_t i = 0; i < ov_top_global.size(); ++i)
            std::cerr << "  " << std::setw(4) << (i+1)
                      << "     0x" << std::hex << std::setw(4) << std::setfill('0') << ov_top_global[i].first
                      << std::dec << std::setfill(' ')
                      << "     " << ov_top_global[i].second << "\n";
    }

    return { total_inserted.load(), total_written.load() };
}
