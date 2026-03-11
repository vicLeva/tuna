#pragma once

// KMC-style partition step — Phase 1 alternative.
//
// Key differences from the cuttlefish (hash-based) and kmtricks (repart) approaches:
//
//   1. Minimizer ordering: uses the raw 2-bit-encoded m-mer value (lex order)
//      rather than a hash, with a validity filter that rejects low-complexity
//      patterns (homopolymers, TGT, TG* suffixes, AA runs).  This is KMC's
//      CMmer::norm[] trick: norm[i] = min(valid(i), valid(rc(i))), where
//      "valid" maps an invalid m-mer to a sentinel value SPECIAL = 4^m.
//
//   2. Window tracking: lazy O(k) rescan only when the current minimizer falls
//      out of the k-mer window, O(1) otherwise.  Mirrors CSplitter in KMC.
//
//   3. Table size: 4^m entries (e.g. m=9 → 262 144 entries, 1 MB) — directly
//      indexed by the raw m-mer integer, unlike the repart table which uses a
//      hash-truncation key over 2^(2*REPART_BITS) entries.
//
//   4. Load-balanced table: a pre-scan (Phase 0) counts the actual number of
//      k-mers whose window-minimizer resolves to each norm value; the LPT
//      greedy assigns norm values to partitions to equalise k-mer load.
//
// Public surface:
//   partition_kmers_kmc<k,l>(cfg, buckets)   — parallel harness
//
// Modifiable bricks:
//   KmcNormTable        — the validity filter and canonicalisation rule.
//   KmcPartMap          — the norm-value → partition mapping (uniform / LPT).
//   extract_superkmers_kmc<k> — pure sequence logic, no I/O, no threads.
//   scan_kmc_sig_counts<k>   — pre-scan counting k-mers per norm value.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "fast_fasta.hpp"

#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <atomic>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <list>
#include <cstdint>


// ─── KMC base encoding ────────────────────────────────────────────────────────
// A=0, C=1, G=2, T=3 — the encoding used by CMmer in KMC.
// Any other character (N, gap, …) → KMC_INVALID.

static constexpr uint32_t KMC_INVALID = 4u;

inline uint32_t kmc_enc(char c)
{
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return KMC_INVALID;
    }
}


// ─── KMC norm table ───────────────────────────────────────────────────────────
//
// Pre-computed lookup: norm[i] = canonical signature for m-mer encoded as i.
//   norm[i] = min( is_allowed(i)     ? i     : SPECIAL,
//                  is_allowed(rc(i)) ? rc(i) : SPECIAL )
//
// SPECIAL = 4^m (one past the last valid encoding).
// Invalid m-mers (both strands rejected) map to SPECIAL; they are assigned to
// the designated "overflow" partition by KmcPartMap.
//
// Validity filter (ported verbatim from KMC's CMmer::is_allowed):
//   • TTT suffix, TGT suffix, TG* suffix rejected
//   • AA di-nucleotide anywhere inside rejected
//   • AAA / ACA / *AA prefix rejected
// This avoids assigning large k-mer populations to a single partition.

struct KmcNormTable
{
    uint32_t m;       // signature length (KMC default: 9; max supported: 11)
    uint32_t special; // = 4^m; sentinel for "no valid signature"
    uint32_t mask;    // = special - 1; used to keep m-mer in [0, special)

    std::vector<uint32_t> norm; // size = special + 1

    explicit KmcNormTable(uint32_t sig_len) : m(sig_len)
    {
        assert(m >= 5 && m <= 11 && "KMC signature length must be in [5, 11]");
        special = 1u << (2 * m);
        mask    = special - 1u;
        norm.resize(special + 1u);

        for (uint32_t i = 0; i < special; ++i) {
            const uint32_t rc   = kmc_rc(i);
            const uint32_t fval = kmc_is_allowed(i)  ? i  : special;
            const uint32_t rval = kmc_is_allowed(rc) ? rc : special;
            norm[i] = std::min(fval, rval);
        }
        norm[special] = special; // SPECIAL maps to itself
    }

    uint32_t operator[](uint32_t mmer) const { return norm[mmer]; }

    // Validity filter, ported verbatim from KMC mmer.h CMmer::is_allowed().
    // Public so KmcPartMap::from_counts can iterate over all allowed raw values.
    bool kmc_is_allowed(uint32_t mmer) const
    {
        if ((mmer & 0x3fu) == 0x3fu) return false; // TTT suffix
        if ((mmer & 0x3fu) == 0x3bu) return false; // TGT suffix
        if ((mmer & 0x3cu) == 0x3cu) return false; // TG* suffix

        uint32_t tmp = mmer;
        for (uint32_t j = 0; j < m - 3; ++j) {
            if ((tmp & 0xfu) == 0u) return false;  // AA inside
            tmp >>= 2;
        }
        if (tmp        == 0u)    return false; // AAA prefix
        if (tmp        == 0x04u) return false; // ACA prefix
        if ((tmp & 0xfu) == 0u) return false;  // *AA prefix

        return true;
    }

private:

    // Reverse complement of a 2-bit-encoded m-mer.
    // Complement: 3-x swaps A↔T, C↔G.  Then reverse the order.
    uint32_t kmc_rc(uint32_t mmer) const
    {
        uint32_t rev   = 0;
        uint32_t shift = 2 * (m - 1);
        for (uint32_t i = 0; i < m; ++i) {
            rev  += (3u - (mmer & 3u)) << shift;
            mmer >>= 2;
            shift -= 2;
        }
        return rev;
    }
};


// ─── Raw m-mer → partition map ────────────────────────────────────────────────
//
// Maps each raw m-mer (0 … SPECIAL) to a partition ID.
// SPECIAL always maps to the last partition (overflow bin).
//
// Two builders:
//   uniform(nt, n)            — norm_val % n (no pre-scan)
//   from_counts(nt, cnts, n)  — exact port of KMC's CSignatureMapper::Init
//
// NOTE: from_counts fills map[raw_mmer] directly (like KMC's signature_map).
// Non-canonical raw values are included in the greedy but are never looked up
// (the partition code always passes a norm/canonical value).

struct KmcPartMap
{
    std::vector<uint32_t> map; // size = special + 1; map[raw_mmer] → partition

    uint32_t operator[](uint32_t raw_mmer) const { return map[raw_mmer]; }

    static KmcPartMap uniform(const KmcNormTable& nt, uint32_t num_parts)
    {
        KmcPartMap pm;
        pm.map.resize(nt.special + 1u);
        for (uint32_t i = 0; i < nt.special; ++i)
            pm.map[i] = (nt.norm[i] == nt.special) ? (num_parts - 1u)
                                                    : (nt.norm[i] % num_parts);
        pm.map[nt.special] = num_parts - 1u; // overflow bin
        return pm;
    }

    // Exact port of KMC's CSignatureMapper::Init (kmc_core/s_mapper.h).
    //
    // Differences from the old LPT approach:
    //   • Includes ALL kmc_is_allowed raw values (not just canonical norm values).
    //     Non-canonical ones have counts[r]=0 from the Phase-0 scan (which always
    //     writes to the canonical/norm index); they receive the +1000 baseline and
    //     sort to the bottom — exactly as in KMC.
    //   • Greedy strategy: fill-up with 1.1×mean tolerance (not LPT min-heap).
    //   • Initial mean = sum / num_parts (KMC divides by n_bins, including the
    //     special/overflow bin, giving a slightly lower mean).
    //   • +1000 baseline per signature (KMC's heuristic to avoid zero-count bins).
    //   • map[] filled by raw value directly (like KMC's signature_map[]).
    static KmcPartMap from_counts(
        const KmcNormTable&          nt,
        const std::vector<uint64_t>& counts, // counts[norm_val] from Phase-0 scan
        uint32_t                     num_parts)
    {
        assert(counts.size() == nt.special + 1u);

        // 1. Collect (raw_val, count) for every allowed m-mer, sorted descending.
        //    counts[r] == 0 for non-canonical r (Phase-0 always writes to norm[r]).
        using Entry = std::pair<uint32_t, uint64_t>; // (raw_val, count)
        std::list<Entry> stats;
        for (uint32_t r = 0; r < nt.special; ++r)
            if (nt.kmc_is_allowed(r))
                stats.push_back({r, counts[r]});

        stats.sort([](const Entry& a, const Entry& b){ return a.second > b.second; });

        // 2. Add +1000 baseline and accumulate sum (KMC does this after sorting).
        double sum = 0.0;
        for (auto& e : stats) {
            e.second += 1000;
            sum += static_cast<double>(e.second);
        }

        // 3. Fill-up greedy — exact replica of KMC's CSignatureMapper::Init loop.
        const uint32_t max_bins   = num_parts - 1u; // bins for valid signatures
        uint32_t       n          = max_bins;        // remaining bins available
        double         mean       = sum / num_parts; // KMC: sum / n_bins (incl. special)
        double         max_bin_sz = 1.1 * mean;
        uint32_t       bin_no     = 0;

        // Default: every raw value → overflow bin (handles filtered / SPECIAL).
        KmcPartMap pm;
        pm.map.assign(nt.special + 1u, num_parts - 1u);

        std::list<Entry> group;

        while (stats.size() > n) {
            Entry& front = stats.front();
            if (static_cast<double>(front.second) > mean) {
                // Heavy: assign its own bin.
                pm.map[front.first] = bin_no++;
                sum       -= static_cast<double>(front.second);
                mean       = sum / (max_bins - bin_no);
                max_bin_sz = 1.1 * mean;
                stats.pop_front();
                --n;
            } else {
                // Fill-up: pack as many entries as fit under max_bin_sz.
                group.clear();
                double tmp_sum = 0.0;
                for (auto it = stats.begin(); it != stats.end(); ) {
                    // Early-exit: if even the smallest entry overflows, stop.
                    if (tmp_sum + static_cast<double>(stats.back().second) >= max_bin_sz)
                        break;
                    if (tmp_sum + static_cast<double>(it->second) < max_bin_sz) {
                        tmp_sum += static_cast<double>(it->second);
                        group.push_back(*it);
                        it = stats.erase(it);
                    } else {
                        ++it;
                    }
                }
                for (const auto& g : group)
                    pm.map[g.first] = bin_no;
                --n;
                ++bin_no;
                sum       -= tmp_sum;
                mean       = sum / (max_bins - bin_no);
                max_bin_sz = 1.1 * mean;
            }
        }

        // Remaining entries: one per bin (KMC's cleanup loop).
        for (const auto& e : stats)
            pm.map[e.first] = bin_no++;

        // Special signature → last bin (bin_no == num_parts - 1 here).
        pm.map[nt.special] = bin_no;

        return pm;
    }
};


// ─── Superkmer extraction brick (KMC-style) ───────────────────────────────────
//
// Walk one sequence, maintaining the minimum-norm m-mer in each k-mer window
// using KMC's lazy rescan strategy (mirrors CSplitter::ProcessReads):
//   • end_mmer slides one base at a time.
//   • If end_mmer < current_sig: new minimum found → superkmer boundary.
//   • If end_mmer == current_sig: update position (keep rightmost occurrence).
//   • If current_sig has fallen out of the k-mer window: rescan window.
// Superkmers are written to the bucket determined by part_map[current_sig].
//
// No I/O, no locking, no thread state: pure, self-contained, testable.

template <uint16_t k>
void extract_superkmers_kmc(
    const char* const              seq,
    const size_t                   seq_len,
    const KmcNormTable&            nt,
    const KmcPartMap&              pm,
    std::vector<SuperkmerWriter>&  writers,
    uint64_t&                      kmer_count)
{
    const uint32_t m    = nt.m;
    const uint32_t mask = nt.mask;

    if (seq_len < k) return;

    // Build a raw m-mer string (uint32_t) from m consecutive characters.
    const auto build_str = [&](size_t pos) -> uint32_t {
        uint32_t s = 0;
        for (uint32_t j = 0; j < m; ++j)
            s = ((s << 2) | kmc_enc(seq[pos + j])) & mask;
        return s;
    };

    size_t i = 0;
    while (i + k <= seq_len) {

        // Find the next position with m consecutive valid (non-N) bases.
        bool found_n = false;
        for (uint32_t j = 0; j < m; ++j) {
            if (kmc_enc(seq[i + j]) == KMC_INVALID) {
                i += j + 1;
                found_n = true;
                break;
            }
        }
        if (found_n) continue;
        if (i + k > seq_len) break;

        // Initialise with the first m-mer at position i.
        uint32_t cur_str = build_str(i);
        uint32_t cur_val = nt[cur_str];
        uint32_t end_str = cur_str;
        uint32_t end_val = cur_val;

        size_t sig_start = i;  // where current_sig starts in the sequence
        size_t sk_start  = i;  // start of the current superkmer
        size_t len       = m;  // accumulated superkmer length
        size_t pos       = i + m;

        while (pos < seq_len) {
            const uint32_t base = kmc_enc(seq[pos]);

            if (base == KMC_INVALID) {
                // Sequence break — emit and restart from after the N.
                if (len >= k) {
                    writers[pm[cur_val]].append(seq + sk_start, len);
                    kmer_count += len - k + 1;
                }
                len = 0;
                ++pos;
                break;
            }

            // Slide end m-mer one base to the right.
            end_str = ((end_str << 2) | base) & mask;
            end_val = nt[end_str];

            if (end_val < cur_val) {
                // New minimum — superkmer boundary.
                if (len >= k) {
                    writers[pm[cur_val]].append(seq + sk_start, len);
                    kmer_count += len - k + 1;
                    sk_start = pos - (k - 1);
                    len      = k - 1;
                }
                cur_str   = end_str; cur_val = end_val;
                sig_start = pos - m + 1;

            } else if (end_val == cur_val) {
                // Tie — prefer the rightmost occurrence (like KMC).
                cur_str   = end_str; cur_val = end_val;
                sig_start = pos - m + 1;

            } else if (sig_start + k - 1 < pos) {
                // Current minimizer has fallen out of the k-mer window.
                // Emit the current superkmer, then rescan the window.
                if (len >= k) {
                    writers[pm[cur_val]].append(seq + sk_start, len);
                    kmer_count += len - k + 1;
                    sk_start = pos - (k - 1);
                    len      = k - 1;
                }
                ++sig_start;

                // Rescan [sig_start … pos] for the new minimum.
                end_str = build_str(sig_start);
                end_val = nt[end_str];
                cur_str = end_str; cur_val = end_val;
                for (size_t j = sig_start + m; j <= pos; ++j) {
                    end_str = ((end_str << 2) | kmc_enc(seq[j])) & mask;
                    end_val = nt[end_str];
                    if (end_val <= cur_val) {
                        cur_str   = end_str; cur_val = end_val;
                        sig_start = j - m + 1;
                    }
                }
            }

            ++len;
            ++pos;
        }

        // Emit the last superkmer of this run.
        if (len >= k) {
            writers[pm[cur_val]].append(seq + sk_start, len);
            kmer_count += len - k + 1;
        }

        i = pos;
    }
}


// ─── Phase 0: KMC signature frequency scan ────────────────────────────────────
//
// Walks all input sequences using the same KMC window algorithm and counts
// how many k-mers resolve to each norm value.  Returns a vector of size
// special + 1 = 4^m + 1; counts[norm_val] = total k-mers whose window
// minimizer (under the KMC norm ordering) equals norm_val.
//
// This is the EXACT count used by KMC's CSignatureMapper::Init().

template <uint16_t k>
std::vector<uint64_t> scan_kmc_sig_counts(
    const Config&       cfg,
    const KmcNormTable& nt)
{
    const uint32_t m       = nt.m;
    const uint32_t mask    = nt.mask;
    const uint32_t table_sz = nt.special + 1u;

    const size_t n_files   = cfg.input_files.size();
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_files);

    std::vector<std::vector<uint64_t>> thr_counts(
        n_threads, std::vector<uint64_t>(table_sz, 0));

    auto worker = [&](size_t tid) {
        auto& local = thr_counts[tid];

        const auto build_str = [&](const char* s) -> uint32_t {
            uint32_t v = 0;
            for (uint32_t j = 0; j < m; ++j)
                v = ((v << 2) | kmc_enc(s[j])) & mask;
            return v;
        };

        SeqReader parser; // one reader per thread, reused across files
        for (size_t fi = tid; fi < n_files; fi += n_threads) {
            if (!parser.load(cfg.input_files[fi])) continue;

            while (parser.read_next_seq()) {
                const char*  seq     = parser.seq();
                const size_t seq_len = parser.seq_len();
                if (seq_len < k) continue;

                size_t i = 0;
                while (i + k <= seq_len) {
                    bool found_n = false;
                    for (uint32_t j = 0; j < m; ++j) {
                        if (kmc_enc(seq[i + j]) == KMC_INVALID) {
                            i += j + 1;
                            found_n = true;
                            break;
                        }
                    }
                    if (found_n) continue;
                    if (i + k > seq_len) break;

                    uint32_t cur_str = build_str(seq + i);
                    uint32_t cur_val = nt[cur_str];
                    uint32_t end_str = cur_str;
                    uint32_t end_val = cur_val;

                    size_t sig_start = i;
                    size_t len       = m;
                    size_t pos       = i + m;

                    while (pos < seq_len) {
                        const uint32_t base = kmc_enc(seq[pos]);

                        if (base == KMC_INVALID) {
                            if (len >= k)
                                local[cur_val] += len - k + 1;
                            len = 0;
                            ++pos;
                            break;
                        }

                        end_str = ((end_str << 2) | base) & mask;
                        end_val = nt[end_str];

                        if (end_val < cur_val) {
                            if (len >= k) {
                                local[cur_val] += len - k + 1;
                                len = k - 1;
                            }
                            cur_str = end_str; cur_val = end_val;
                            sig_start = pos - m + 1;

                        } else if (end_val == cur_val) {
                            cur_str = end_str; cur_val = end_val;
                            sig_start = pos - m + 1;

                        } else if (sig_start + k - 1 < pos) {
                            if (len >= k) {
                                local[cur_val] += len - k + 1;
                                len = k - 1;
                            }
                            ++sig_start;
                            end_str = build_str(seq + sig_start);
                            end_val = nt[end_str];
                            cur_str = end_str; cur_val = end_val;
                            for (size_t j = sig_start + m; j <= pos; ++j) {
                                end_str = ((end_str << 2) | kmc_enc(seq[j])) & mask;
                                end_val = nt[end_str];
                                if (end_val <= cur_val) {
                                    cur_str = end_str; cur_val = end_val;
                                    sig_start = j - m + 1;
                                }
                            }
                        }

                        ++len;
                        ++pos;
                    }

                    if (len >= k)
                        local[cur_val] += len - k + 1;

                    i = pos;
                }
            }

        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    std::vector<uint64_t> counts(table_sz, 0);
    for (size_t t = 0; t < n_threads; ++t)
        for (uint32_t i = 0; i < table_sz; ++i)
            counts[i] += thr_counts[t][i];

    return counts;
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// Spawns min(num_threads, n_files) threads, round-robin file assignment.
// Each thread runs extract_superkmers_kmc per sequence, flushing per-bucket
// SuperkmerWriters to the shared ofstreams under per-bucket mutexes.
//
// The template parameter l is carried for interface consistency with the other
// partition_kmers variants; the KMC signature length comes from cfg.kmc_sig_len.

template <uint16_t k, uint16_t /* l */>
PartitionStats partition_kmers_kmc(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    const KmcNormTable&         nt,
    const KmcPartMap&           pm)
{
    const size_t n_files   = cfg.input_files.size();
    const size_t n_threads = std::min(static_cast<size_t>(cfg.num_threads), n_files);

    std::vector<std::mutex> bucket_mutexes(cfg.num_partitions);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0};

    auto worker = [&](size_t tid) {
        std::vector<SuperkmerWriter> writers(cfg.num_partitions);
        uint64_t local_seqs = 0, local_kmers = 0;

        SeqReader parser; // one reader per thread, reused across files
        for (size_t fi = tid; fi < n_files; fi += n_threads) {
            if (!parser.load(cfg.input_files[fi])) continue;

            while (parser.read_next_seq()) {
                ++local_seqs;

                extract_superkmers_kmc<k>(
                    parser.seq(), parser.seq_len(),
                    nt, pm, writers, local_kmers);

                for (size_t p = 0; p < cfg.num_partitions; ++p)
                    if (writers[p].needs_flush())
                        writers[p].flush_to(buckets[p], bucket_mutexes[p]);
            }

            parser.close();
        }

        for (size_t p = 0; p < cfg.num_partitions; ++p)
            writers[p].flush_to(buckets[p], bucket_mutexes[p]);

        total_seqs .fetch_add(local_seqs,  std::memory_order_relaxed);
        total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    return { total_seqs.load(), total_kmers.load() };
}
