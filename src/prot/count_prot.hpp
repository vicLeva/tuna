#pragma once

// Phase 2+3 — superkmer partitions → counted protein k-mers → TSV output.
//
// Hash table: flat open-addressing (linear probe), key=uint64_t, value=uint32_t.
// Protein k-mers (k ≤ 12) pack into uint64_t at 5 bits/AA.
// Sentinel: UINT64_MAX (valid k-mers use at most 60 bits for k=12).
//
// Prefetch: 1-ahead prefetch of the next superkmer's first k-mer slot hides
// ~40 ns LLC miss behind the current superkmer's counting loop.

#include "config_prot.hpp"
#include "superkmer_prot.hpp"
#include "aa_encode.hpp"
#include "aa_hash.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// xxHash for hash-table bucket routing.
#include "xxHash/xxhash.h"

namespace prot {

// ── Protein k-mer hash table ──────────────────────────────────────────────────
//
// Linear-probe open-addressing with power-of-2 capacity.
// Empty slot sentinel: EMPTY_KEY = UINT64_MAX.
// Valid protein k-mers for k ≤ 12 never reach UINT64_MAX (max = (1<<60)-1).

struct ProtKmerTable {
    static constexpr uint64_t EMPTY_KEY = ~uint64_t(0);

    struct alignas(8) Entry {
        uint64_t key   = EMPTY_KEY;
        uint32_t count = 0;
        uint32_t _pad  = 0;
    };

    Entry*  table_ = nullptr;
    size_t  mask_  = 0;
    size_t  size_  = 0;

    ProtKmerTable() = default;
    explicit ProtKmerTable(size_t min_cap) {
        size_t cap = 1;
        while (cap < min_cap) cap <<= 1;
        // Over-provision for load factor: target 70% max → allocate cap / 0.7.
        cap = std::max(cap, size_t(64));
        table_ = new Entry[cap]();
        mask_  = cap - 1;
        size_  = 0;
    }

    ProtKmerTable(const ProtKmerTable&)            = delete;
    ProtKmerTable& operator=(const ProtKmerTable&) = delete;

    ProtKmerTable(ProtKmerTable&& o) noexcept
        : table_(o.table_), mask_(o.mask_), size_(o.size_)
    { o.table_ = nullptr; o.mask_ = 0; o.size_ = 0; }

    ~ProtKmerTable() { delete[] table_; }

    void clear() {
        std::fill(table_, table_ + mask_ + 1, Entry{EMPTY_KEY, 0, 0});
        size_ = 0;
    }

    size_t capacity() const noexcept { return mask_ + 1; }

    __attribute__((always_inline))
    void upsert(uint64_t key) noexcept {
        size_t h = XXH3_64bits(&key, 8) & mask_;
        while (true) {
            Entry& e = table_[h];
            if (__builtin_expect(e.key == key, 1)) { ++e.count; return; }
            if (e.key == EMPTY_KEY) { e.key = key; e.count = 1; ++size_; return; }
            h = (h + 1) & mask_;
        }
    }

    // Prefetch the slot for key without reading/writing it.
    __attribute__((always_inline))
    void prefetch(uint64_t key) const noexcept {
        const size_t h = XXH3_64bits(&key, 8) & mask_;
        __builtin_prefetch(&table_[h], 1, 3);
    }

    // Iterate over all entries.
    template <typename F>
    void for_each(F&& f) const {
        for (size_t i = 0; i <= mask_; ++i)
            if (table_[i].key != EMPTY_KEY)
                f(table_[i].key, table_[i].count);
    }
};


// ── Unpack helpers ────────────────────────────────────────────────────────────

// Decode the first k amino acids from 5-bit packed data into an array.
template <uint16_t k>
inline void unpack_k_aas(const uint8_t* packed, uint8_t* out) noexcept {
    for (uint32_t i = 0; i < k; ++i)
        out[i] = read5(packed, i * 5);
}

// Build a k-mer uint64_t from k consecutive encoded amino acids.
template <uint16_t k>
inline uint64_t build_kmer(const uint8_t* enc) noexcept {
    uint64_t kmer = 0;
    for (uint16_t i = 0; i < k; ++i)
        kmer = (kmer << 5) | enc[i];
    return kmer;
}

// Compute the first k-mer from 5-bit packed data (first k amino acids).
template <uint16_t k>
inline uint64_t first_kmer_from_packed(const uint8_t* packed) noexcept {
    uint64_t kmer = 0;
    for (uint16_t i = 0; i < k; ++i)
        kmer = (kmer << 5) | read5(packed, i * 5);
    return kmer;
}

static constexpr uint64_t KMER_MASK_FN(uint16_t k) noexcept {
    return (k >= 12) ? uint64_t(0x0FFFFFFFFFFFFFFF) : ((uint64_t(1) << (5 * k)) - 1);
}


// ── Counting brick ────────────────────────────────────────────────────────────
//
// Drains a superkmer reader into a ProtKmerTable.
// 1-ahead prefetch hides the LLC miss for the next superkmer's first k-mer slot.

template <uint16_t k, uint16_t m, typename Reader = ProtSuperKmerReader<k, m>>
uint64_t count_prot_partition(Reader& reader, ProtKmerTable& table) {
    [[maybe_unused]] static constexpr auto NO_MIN = psk_no_min<k, m>;
    static constexpr uint64_t kmer_mask = KMER_MASK_FN(k);

    uint64_t inserted = 0;
    uint8_t  dec_buf[512];  // decoded AAs for current superkmer (max superkmer = 2k-m)

    if (!reader.next()) return inserted;

    const uint8_t* cur_packed  = reader.packed_data();
    size_t         cur_len     = reader.size();

    // Prefetch loop: prime the pump, then 1-ahead prefetch.
    while (cur_len >= k) {
        // Prefetch next superkmer's first k-mer slot.
        const uint8_t* next_packed = nullptr;
        size_t         next_len    = 0;
        bool           has_next    = false;

        if (reader.next()) {
            next_packed = reader.packed_data();
            next_len    = reader.size();
            has_next    = true;
            if (next_len >= k) {
                const uint64_t nk = first_kmer_from_packed<k>(next_packed);
                table.prefetch(nk);
            }
        }

        // Decode current superkmer's amino acids.
        const size_t n = cur_len;
        for (size_t i = 0; i < n; ++i)
            dec_buf[i] = read5(cur_packed, static_cast<uint32_t>(i * 5));

        // Build initial k-mer and slide.
        uint64_t kmer = build_kmer<k>(dec_buf);
        table.upsert(kmer);
        ++inserted;

        for (size_t pos = k; pos < n; ++pos) {
            kmer = ((kmer << 5) | dec_buf[pos]) & kmer_mask;
            table.upsert(kmer);
            ++inserted;
        }

        cur_packed = next_packed;
        cur_len    = next_len;
        if (!has_next) break;
    }

    return inserted;
}


// ── Output brick ─────────────────────────────────────────────────────────────

template <uint16_t k>
void write_prot_counts(const ProtKmerTable& table,
                       std::ostream& out,
                       std::mutex& out_mtx) {
    // Buffer output to avoid per-entry lock contention.
    std::string buf;
    buf.reserve(256 << 10);
    char kmer_str[k + 2];

    table.for_each([&](uint64_t kmer, uint32_t count) {
        decode_pkmer<k>(kmer, kmer_str);
        buf += kmer_str;
        buf += '\t';
        buf += std::to_string(count);
        buf += '\n';
        if (buf.size() > (128u << 10)) {
            std::lock_guard<std::mutex> g(out_mtx);
            out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
            buf.clear();
        }
    });
    if (!buf.empty()) {
        std::lock_guard<std::mutex> g(out_mtx);
        out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
    }
}


// ── Parallel harness: in-memory ───────────────────────────────────────────────

template <uint16_t k, uint16_t m>
void count_and_write_prot_mem(
    const ProtConfig&         cfg,
    const std::vector<std::string>& part_bufs,
    std::ostream&             out)
{
    const size_t n_parts   = part_bufs.size();
    const size_t n_threads = cfg.num_threads;

    std::atomic<size_t> part_idx{0};
    std::mutex          out_mtx;
    std::exception_ptr  worker_ex;
    std::mutex          ex_mtx;

    // init_sz heuristic: 2 × estimated_kmers_per_partition / coverage
    // For proteins: use partition buffer size / bytes_per_aa_avg as estimate.
    const uint64_t total_buf = [&]() {
        uint64_t s = 0;
        for (auto& b : part_bufs) s += b.size();
        return s;
    }();
    // Each entry in the buffer is ~HDR_BYTES + 5*sk_len/8 bytes; rough estimate:
    // unique k-mers ≈ total_buf * 8 / 5 / k (very rough).
    const size_t init_sz_per_part = [&]() {
        if (n_parts == 0) return size_t(1 << 14);
        const size_t est = static_cast<size_t>(total_buf * 8 / 5 / cfg.k / n_parts * 2);
        return std::clamp(est, size_t(1 << 12), size_t(1 << 22));
    }();

    auto worker = [&]() {
        try {
            for (;;) {
                const size_t pi = part_idx.fetch_add(1);
                if (pi >= n_parts) break;
                if (part_bufs[pi].empty()) continue;

                ProtMemoryReader<k, m> reader(part_bufs[pi]);
                ProtKmerTable table(init_sz_per_part);
                count_prot_partition<k, m, ProtMemoryReader<k, m>>(reader, table);
                write_prot_counts<k>(table, out, out_mtx);
            }
        } catch (...) {
            std::lock_guard<std::mutex> g(ex_mtx);
            if (!worker_ex) worker_ex = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t i = 0; i < n_threads; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    if (worker_ex) std::rethrow_exception(worker_ex);
}


// ── Parallel harness: disk ────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
void count_and_write_prot_disk(
    const ProtConfig& cfg,
    const std::string& work_dir,
    size_t             n_parts,
    std::ostream&      out)
{
    const size_t n_threads = cfg.num_threads;
    std::atomic<size_t> part_idx{0};
    std::mutex          out_mtx;
    std::exception_ptr  worker_ex;
    std::mutex          ex_mtx;

    constexpr size_t DEFAULT_INIT_SZ = 1 << 16;

    auto worker = [&]() {
        try {
            for (;;) {
                const size_t pi = part_idx.fetch_add(1);
                if (pi >= n_parts) break;
                const std::string path = prot_partition_path(work_dir, pi);
                ProtSuperKmerReader<k, m> reader(path);
                if (!reader.ok()) continue;
                ProtKmerTable table(DEFAULT_INIT_SZ);
                count_prot_partition<k, m>(reader, table);
                write_prot_counts<k>(table, out, out_mtx);
            }
        } catch (...) {
            std::lock_guard<std::mutex> g(ex_mtx);
            if (!worker_ex) worker_ex = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t i = 0; i < n_threads; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    if (worker_ex) std::rethrow_exception(worker_ex);
}

} // namespace prot
