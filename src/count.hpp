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
#include <unordered_map>

#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"


// ─── Debug stats ──────────────────────────────────────────────────────────────
//
// Populated by count_partition when dbg != nullptr.
// coverage_hist: maps (k-mers sharing one minimizer) → (number of minimizers
// with that coverage).  Aggregated globally in count_and_write and written to
// <work_dir>/debug_min_coverage.csv.

struct PartitionDebugInfo {
    uint64_t n_inserted  = 0;   // total k-mer insertions (with multiplicity)
    uint64_t n_overflow  = 0;   // insertions that went to overflow table
    uint64_t table_cap   = 0;   // final flat-table capacity (k-mer slots)
    uint64_t n_buckets   = 0;   // table_cap / B

    // coverage → count_of_minimizers_with_that_coverage
    // Only populated for superkmers with min_pos != 0xFF.
    std::unordered_map<uint32_t, uint64_t> coverage_hist;
};


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
    typename kache_hash::Streaming_Kmer_Hash_Table<k, false, uint32_t, l>::Token& token,
    PartitionDebugInfo* dbg = nullptr)
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

    // Per-minimizer k-mer count (only allocated when debug is requested).
    // Maps minimizer_hash → total k-mers sharing that minimizer in this partition.
    std::unordered_map<uint64_t, uint32_t> min_kmer_count;

    while (reader.next()) {
        const uint8_t* packed   = reader.packed_data();
        const size_t   len      = reader.size();
        const uint8_t  min_pos  = reader.min_pos();
        if (len < k) continue;

        const uint32_t sk_kmers = static_cast<uint32_t>(len - k + 1);

        // min_pos == 0xFF is a sentinel (KMC mode: multiple ntHash minimizers
        // possible within one superkmer) — fall back to MinimizerWindow::reset().
        // Otherwise use the stored position to compute ntHash in O(l) vs O(k).
        kache_hash::Kmer_Window<k, l> win;
        if (min_pos != 0xFF) {
            const uint64_t mh = lmer_nt_hash(packed, min_pos);
            win.init_packed_with_hash(packed, mh);
            if (dbg) min_kmer_count[mh] += sk_kmers;
        } else {
            win.init_packed(packed);
        }

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

    if (dbg) {
        // Convert per-minimizer k-mer counts to a coverage histogram.
        for (auto& [mh, cnt] : min_kmer_count)
            dbg->coverage_hist[cnt]++;
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
// init_sz scales with partition count: 128M / n_parts, clamped to [256K, 4M].
// Two competing constraints drive the formula:
//   (1) Per-bucket occupancy ≤ ~50%: secondary bucket probing triggers when the
//       primary bucket (B=32 slots) is full.  Secondary probe = 2 XXH3 + 2 cold
//       cache line accesses.  Keeping occupancy low avoids this path.
//   (2) Metadata fits in L3: each bucket = 64 B metadata (cs[32] + min_coord[32]).
//       With 4 active threads, total metadata budget = L3/n_threads ≈ 7.5 MB.
// Formula target: ~37% per-bucket occupancy at expected unique k-mers / partition.
//   n=32:   4M → 128K buckets →  8 MB meta (slightly above L3, same as original)
//   n=64:   2M →  64K buckets →  4 MB meta (fits in L3)
//   n=128:  1M →  32K buckets →  2 MB meta
//   n=256: 512K→  16K buckets →  1 MB meta

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

    // Debug statistics (only used when cfg.debug_stats).
    std::mutex                                     dbg_mutex;
    std::unordered_map<uint32_t, uint64_t>         global_coverage_hist;
    std::vector<PartitionDebugInfo>                part_infos;
    if (cfg.debug_stats) part_infos.resize(n_parts);

    auto worker = [&](size_t tid) {
        typename table_t::Token token;
        std::string chunk;

        for (size_t p = tid; p < n_parts; p += n_threads) {
            const std::string part_path = cfg.work_dir + cfg.partition_prefix()
                                          + "_" + std::to_string(p) + ".superkmers";
            SuperkmerReader reader(part_path);
            // Scale initial table size with partition count.
            // Target ~37% per-bucket occupancy to avoid secondary bucket probing
            // (buckets have B=32 slots; secondary probe triggers when primary is full).
            // Also keeps metadata (64 B/bucket) fitting in L3 (30 MB / 4 threads = 7.5 MB).
            // Formula: 128M / n_parts, clamped to [256K, 4M].
            // For n=32:   4M (cap=128K, ~37% occ, 8 MB meta — original, no regression)
            // For n=64:   2M (cap= 64K, ~37% occ, 4 MB meta — fits in L3)
            // For n=128:  1M (cap= 32K, ~37% occ, 2 MB meta)
            // For n=256: 512K(cap= 16K, ~37% occ, 1 MB meta)
            const size_t init_sz = std::clamp(
                size_t(1u << 27) / n_parts,
                size_t(1u << 18),   // min 256K
                size_t(1u << 22));  // max 4M
            table_t table(init_sz, 1);

            PartitionDebugInfo* dbg = cfg.debug_stats ? &part_infos[p] : nullptr;
            const uint64_t ins = count_partition<k, l>(reader, table, token, dbg);
            total_inserted.fetch_add(ins, std::memory_order_relaxed);

            if (dbg) {
                dbg->n_inserted = ins;
                dbg->n_overflow = table.overflow_insert_count();
                dbg->table_cap  = table.capacity();
                dbg->n_buckets  = table.bucket_count();
            }

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

            // Merge debug coverage histogram under lock.
            if (dbg && !dbg->coverage_hist.empty()) {
                std::lock_guard<std::mutex> lg(dbg_mutex);
                for (auto& [cov, cnt] : dbg->coverage_hist)
                    global_coverage_hist[cov] += cnt;
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

        std::cerr << "[overflow] top minimizer bins (top " << table_t::OV_HIST_BITS_PUBLIC
                  << " bits of ntHash canonical):\n";
        std::cerr << "  rank     bin_id (hex)     overflow_kmers\n";
        for(std::size_t i = 0; i < ov_top_global.size(); ++i)
            std::cerr << "  " << std::setw(4) << (i+1)
                      << "     0x" << std::hex << std::setw(4) << std::setfill('0') << ov_top_global[i].first
                      << std::dec << std::setfill(' ')
                      << "     " << ov_top_global[i].second << "\n";
    }

    // ── Debug output ──────────────────────────────────────────────────────────
    if (cfg.debug_stats) {
        // Per-partition table summary.
        std::cerr << "\n[debug] per-partition table stats (first 20):\n";
        std::cerr << "  part   n_inserted   n_overflow   ov%      table_cap   n_buckets   avg_k/bucket\n";
        for (size_t p = 0; p < std::min(n_parts, size_t(20)); ++p) {
            const auto& d = part_infos[p];
            const double ov_pct = d.n_inserted ? 100.0 * d.n_overflow / d.n_inserted : 0.0;
            const double avg_kb = d.n_buckets  ? static_cast<double>(d.n_inserted) / d.n_buckets : 0.0;
            std::cerr << "  " << std::setw(5) << p
                      << "  " << std::setw(11) << d.n_inserted
                      << "  " << std::setw(11) << d.n_overflow
                      << "  " << std::fixed << std::setprecision(1) << std::setw(6) << ov_pct << "%"
                      << "  " << std::setw(10) << d.table_cap
                      << "  " << std::setw(10) << d.n_buckets
                      << "  " << std::setprecision(1) << std::setw(12) << avg_kb << "\n";
        }
        if (n_parts > 20) std::cerr << "  ... (" << (n_parts - 20) << " more partitions)\n";

        // Global minimizer coverage histogram summary.
        if (!global_coverage_hist.empty()) {
            uint64_t total_minimizers = 0, total_kmers_covered = 0;
            uint32_t max_cov = 0;
            for (auto& [cov, cnt] : global_coverage_hist) {
                total_minimizers   += cnt;
                total_kmers_covered += static_cast<uint64_t>(cov) * cnt;
                if (cov > max_cov) max_cov = cov;
            }
            const double avg_cov = total_minimizers ? static_cast<double>(total_kmers_covered) / total_minimizers : 0.0;

            std::cerr << "\n[debug] minimizer coverage (k-mers sharing one minimizer):\n";
            std::cerr << "  unique minimizers tracked : " << total_minimizers << "\n";
            std::cerr << "  total k-mers covered      : " << total_kmers_covered << "\n";
            std::cerr << "  avg k-mers / minimizer    : " << std::fixed << std::setprecision(2) << avg_cov << "\n";
            std::cerr << "  max k-mers / minimizer    : " << max_cov << "\n";

            // Write full histogram to CSV.
            const std::string csv_path = cfg.work_dir + "debug_min_coverage.csv";
            std::ofstream csv(csv_path);
            if (csv) {
                csv << "coverage,n_minimizers,total_kmers\n";
                // Sort by coverage for readability.
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

    return { total_inserted.load(), total_written.load() };
}
