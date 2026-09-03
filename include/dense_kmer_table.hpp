#pragma once
// Dense_Kmer_Hash_Table — phase-2 counting table backed by ankerl::unordered_dense.
//
// This is the phase-2 counting table. It implements the subset of
// kache_hash::Streaming_Kmer_Hash_Table that src/count.hpp uses, which is what
// let the two be benchmarked against each other; unordered_dense won on every
// dataset in the paper suite and became the default. The quotienting table is
// preserved on branch feat/kmer-routing.
//
// ---------------------------------------------------------------------------
// The trade-off this exists to measure
// ---------------------------------------------------------------------------
// The kache table does not store k-mers. Each slot holds one checksum byte and
// one minimizer-coordinate byte, and the key is recovered from the bucket
// address — quotienting. That is roughly 6 bytes per k-mer including the
// count.
//
// unordered_dense must store the key in full: a dense array of
// pair<Kmer<k>, uint32_t> (16 bytes once padded, for k <= 32) plus a bucket
// array of 8 bytes per slot at its 0.8 max load factor, so ~26 bytes per
// k-mer. Expect roughly 3-4x the phase-2 resident memory.
//
// What it buys, and the reason to measure rather than assume: lookups touch a
// contiguous value array with no overflow table, no resize log, no minimizer
// bookkeeping, and no per-thread token state.
//
// ---------------------------------------------------------------------------
// Notes
// ---------------------------------------------------------------------------
// - Kmer<k> is trivially copyable, has operator==, and hashes through to_u64(),
//   so it keys the map directly. Nothing in the output path changes:
//   for_each still yields pair<Kmer<k>, uint32_t> and get_label() still works.
// - prefetch() is a no-op. unordered_dense does not expose its bucket array,
//   so the software-pipelined prefetch that hides the LLC miss in the kache
//   table has no equivalent here. This is a real handicap on large tables and
//   is the first thing to revisit if the measurements are close.
// - Single-threaded per partition only, like the kache table with mt_=false.

#include "kache-hash/Kmer.hpp"
#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"   // ResizeEvent, Kmer_Window

#include <ankerl/unordered_dense.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace tuna {

// to_u64() is already a well-mixed 64-bit hash, so tell unordered_dense not to
// mix it a second time.
template <uint16_t k>
struct Kmer_Hash
{
    using is_avalanching = void;
    uint64_t operator()(const kache_hash::Kmer<k>& key) const noexcept { return key.to_u64(); }
};

template <uint16_t k, bool mt_, typename T_ = uint32_t, uint16_t l = 19, bool canonical_ = true>
class Dense_Kmer_Hash_Table
{
    static_assert(!mt_, "Dense_Kmer_Hash_Table is single-threaded; one table per partition.");

public:
    using key_t = kache_hash::Kmer<k>;
    using val_t = T_;
    using map_t = ankerl::unordered_dense::map<key_t, val_t, Kmer_Hash<k>>;

    // The kache table hands each thread a Token carrying its slot id and
    // resize checkpoint. Nothing here is shared, so the token is empty and
    // exists only to keep the call sites identical.
    struct Token {};

    // count.hpp prints this when reporting overflow bins. There is no overflow
    // table here, so the value is only a label for an always-empty report.
    static constexpr std::size_t OV_HIST_BITS_PUBLIC = 8;

private:
    map_t map_;

public:
    // `resize_worker_c` and `lf` are accepted and ignored: unordered_dense
    // manages its own growth. Reserving up front still matters, since it is
    // what avoids rehashing a partition mid-count.
    explicit Dense_Kmer_Hash_Table(std::size_t max_sz = 1 << 22,
                                   uint64_t resize_worker_c = 1,
                                   double lf = 0.8)
    {
        (void)resize_worker_c; (void)lf;
        map_.reserve(max_sz);
    }

    Dense_Kmer_Hash_Table(const Dense_Kmer_Hash_Table&) = delete;
    Dense_Kmer_Hash_Table& operator=(const Dense_Kmer_Hash_Table&) = delete;

    // Adds the window's current k-mer. `f` is a transform, not a mutator:
    // count.hpp passes `[](uint32_t v){ return v + 1; }`, so an existing count
    // becomes f(count) and a fresh key is stored as `val`. Returns the count
    // before the update, 0 for a first insertion, matching the kache table.
    template <typename F, typename V>
    val_t upsert(const kache_hash::Kmer_Window<k, l>& w, const F& f, const V& val, const Token&)
    {
        const auto& vtx = w.vertex();
        const key_t key = canonical_ ? vtx.canonical() : vtx.kmer();
        const auto [it, inserted] = map_.try_emplace(key, static_cast<val_t>(val));
        if (inserted) return val_t(0);
        const val_t old = it->second;
        it->second = static_cast<val_t>(f(old));
        return old;
    }

    // No equivalent is possible: see the header comment.
    void prefetch(const kache_hash::Kmer_Window<k, l>&) const noexcept {}
    void prefetch_packed(const uint8_t*, uint16_t) const noexcept {}

    template <typename F>
    void for_each(F&& f) const { for (const auto& e : map_) f(e); }

    std::size_t size() const { return map_.size(); }
    std::size_t capacity() const { return map_.bucket_count(); }
    std::size_t bucket_count() const { return map_.bucket_count(); }

    // There is no overflow table and no resize log; the debug reporting in
    // count.hpp reads these unconditionally, so they answer "nothing here".
    uint64_t overflow_insert_count() const { return 0; }
    std::vector<kache_hash::ResizeEvent> resize_log() const { return {}; }
    std::vector<std::pair<uint32_t, uint64_t>> overflow_top_minimizers(std::size_t) const { return {}; }
};

} // namespace tuna
