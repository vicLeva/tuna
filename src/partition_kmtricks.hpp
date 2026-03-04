#pragma once

// kmtricks-style minimizer-to-partition repartition table.
//
// A RepartitionTable maps the top (2*REPART_BITS) bits of a 64-bit
// minimizer hash to a partition ID.  It can be built uniformly
// (equivalent to hash % num_partitions) or load-balanced from a per-key
// frequency scan of the input (Phase 0).
//
// Workflow:
//   // (a) uniform — no pre-scan, replace current approach directly
//   auto rt = RepartitionTable::uniform(cfg.num_partitions);
//
//   // (b) load-balanced — scan first, then build
//   auto counts = scan_minimizer_counts<k>(cfg);
//   auto rt     = RepartitionTable::from_counts(counts, cfg.num_partitions);
//
// Modifiable brick:
//   REPART_BITS        — granularity of the table (4^REPART_BITS entries).
//   RepartitionTable   — the mapping itself; swap the builder to change
//                        the load-balancing algorithm.
//   scan_minimizer_counts — the pre-scan; swap to change what gets counted.

#include "Config.hpp"
#include "fast_fasta.hpp"
#include "DNA_Utility.hpp"
#include "Minimizer_Iterator.hpp"

#include <vector>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <queue>
#include <thread>


// Number of bits used per nucleotide for the table key.
// Table size = 4^REPART_BITS = 2^(2*REPART_BITS) entries.
//   REPART_BITS=10 → 1 048 576 entries, 4 MB per thread during scan.
//   REPART_BITS=12 → 16 777 216 entries, 64 MB per thread (may be heavy).
// Must satisfy: 4^REPART_BITS >= cfg.num_partitions.
static constexpr int REPART_BITS = 10;


// ─── Table data structure ─────────────────────────────────────────────────────

struct RepartitionTable
{
    int                   nb_bits; // bits used; key = hash >> (64 - 2*nb_bits)
    std::vector<uint32_t> table;   // table[key] → partition id

    // Partition lookup: top (2*nb_bits) bits of a 64-bit minimizer hash.
    uint32_t operator()(uint64_t hash) const
    {
        return table[static_cast<size_t>(hash >> (64 - 2 * nb_bits))];
    }

    // Build a uniform table: key % num_parts.
    // Equivalent to the current hash-modulo approach, just via a lookup.
    static RepartitionTable uniform(uint32_t num_parts, int nb_bits = REPART_BITS)
    {
        const uint64_t sz = uint64_t(1) << (2 * nb_bits);
        assert(sz >= num_parts && "REPART_BITS too small for num_partitions");

        RepartitionTable rt;
        rt.nb_bits = nb_bits;
        rt.table.resize(sz);
        for (uint64_t i = 0; i < sz; ++i)
            rt.table[i] = static_cast<uint32_t>(i % num_parts);
        return rt;
    }

    // Build a load-balanced table using LPT (Longest Processing Time) greedy.
    // counts[i] = number of k-mers whose minimizer hash key equals i.
    // Assigns table entries to the least-loaded partition first, aiming to
    // equalise total k-mer count per partition.
    static RepartitionTable from_counts(
        const std::vector<uint64_t>& counts,
        uint32_t                     num_parts,
        int                          nb_bits = REPART_BITS)
    {
        const uint64_t sz = uint64_t(1) << (2 * nb_bits);
        assert(counts.size() == sz);
        assert(sz >= num_parts && "REPART_BITS too small for num_partitions");

        // Sort keys by frequency descending.
        std::vector<uint64_t> order(sz);
        std::iota(order.begin(), order.end(), uint64_t(0));
        std::sort(order.begin(), order.end(),
                  [&](uint64_t a, uint64_t b){ return counts[a] > counts[b]; });

        // Min-heap over (current_load, partition_id).
        using P = std::pair<uint64_t, uint32_t>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> heap;
        for (uint32_t p = 0; p < num_parts; ++p)
            heap.push({uint64_t(0), p});

        RepartitionTable rt;
        rt.nb_bits = nb_bits;
        rt.table.resize(sz);

        for (uint64_t key : order) {
            auto [load, p] = heap.top(); heap.pop();
            rt.table[key]  = p;
            heap.push({load + counts[key], p});
        }
        return rt;
    }
};


// ─── Phase 0: minimizer frequency scan ───────────────────────────────────────
//
// Walk all input sequences with Min_Iterator<k> and count, for every k-mer
// window, how many times each table key (= minimizer hash >> (64 - 2*nb_bits))
// is seen.  Returns a vector of 4^nb_bits counts.
//
// Threading model: min(num_threads, n_files) threads, round-robin file
// assignment, per-thread local uint32_t counts merged to uint64_t at the end.
// Per-thread memory: 4^nb_bits × 4 bytes (≈4 MB for nb_bits=10).

template <uint16_t k>
std::vector<uint64_t> scan_minimizer_counts(const Config& cfg, int nb_bits = REPART_BITS)
{
    const uint64_t table_sz = uint64_t(1) << (2 * nb_bits);
    const int      shift    = 64 - 2 * nb_bits;

    const size_t n_files   = cfg.input_files.size();
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_files);

    // Per-thread uint32_t local counts (4 MB each at nb_bits=10).
    std::vector<std::vector<uint32_t>> thr_counts(
        n_threads, std::vector<uint32_t>(table_sz, 0));

    auto worker = [&](size_t tid) {
        auto& local = thr_counts[tid];
        cuttlefish::Min_Iterator<k> min_it(cfg.l);

        SeqReader parser;
        for (size_t fi = tid; fi < n_files; fi += n_threads) {
            if (!parser.load(cfg.input_files[fi])) continue;

            while (parser.read_next_seq()) {
                const char*  seq     = parser.seq();
                const size_t seq_len = parser.seq_len();
                if (seq_len < k) continue;

                bool in_run = false;

                for (size_t pos = 0; pos < seq_len; ++pos) {
                    const char ch = seq[pos];

                    if (cuttlefish::DNA_Utility::is_placeholder(ch)) {
                        in_run = false;
                        continue;
                    }

                    if (!in_run) {
                        if (pos + k > seq_len) break;

                        bool ok = true;
                        for (size_t t = 1; t < k; ++t)
                            if (cuttlefish::DNA_Utility::is_placeholder(seq[pos + t])) {
                                pos += t; ok = false; break;
                            }
                        if (!ok) continue;

                        min_it.reset(seq + pos);
                        local[static_cast<size_t>(min_it.hash() >> shift)]++;
                        in_run = true;
                        pos   += k - 1;  // for-loop will ++pos to pos+k
                        continue;
                    }

                    min_it.advance(ch);
                    local[static_cast<size_t>(min_it.hash() >> shift)]++;
                }
            }

        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    // Merge per-thread uint32_t counts into a single uint64_t array.
    std::vector<uint64_t> counts(table_sz, 0);
    for (size_t t = 0; t < n_threads; ++t)
        for (uint64_t i = 0; i < table_sz; ++i)
            counts[i] += thr_counts[t][i];

    return counts;
}
