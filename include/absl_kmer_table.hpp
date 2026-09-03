#pragma once
// Absl_Kmer_Hash_Table — phase-2 counting table backed by absl::flat_hash_map.
//
// Third option alongside kache-hash and unordered_dense, selected with:
//
//     cmake -DTUNA_PHASE2_ABSL=ON ...
//
// ---------------------------------------------------------------------------
// Why flat_hash_map and not the others in the library
// ---------------------------------------------------------------------------
// absl offers flat_hash_map, node_hash_map, and the btree maps. For a
// Kmer<k>/uint32_t pair — trivially copyable, 12 bytes before padding — only
// flat_hash_map is a serious candidate:
//   node_hash_map allocates each element separately to give pointer stability
//     that nothing here needs, and pays an indirection on every lookup.
//   btree_map is ordered and comparison-based: O(log n) probes against O(1),
//     worth it only when iteration order matters, which it does not here.
// flat_hash_map is the SwissTable: open addressing, one control byte per slot
// holding a 7-bit hash fragment, probed 16 slots at a time with SSE2. On a
// miss it usually touches one cache line of control bytes and no slot at all.
//
// ---------------------------------------------------------------------------
// Expected footprint
// ---------------------------------------------------------------------------
// One slot is sizeof(pair<Kmer<k>, uint32_t>) = 16 B for k <= 32, plus one
// control byte, at a max load factor of 7/8. So ~19.4 B per stored k-mer at
// full load, against ~26 B for unordered_dense (8 B bucket + 16 B value at
// 0.8) and ~6 B for the kache quotienting table.
//
// Capacity is always a power of two minus one slot, and it grows by doubling,
// so a table sized for n elements can hold up to 2n before rehashing. Reserve
// up front for the same reason as the other two: a mid-count rehash on a large
// partition is the expensive event.
//
// ---------------------------------------------------------------------------
// Notes
// ---------------------------------------------------------------------------
// - Kmer<k>::to_u64() is XXH3, so the hash is already avalanched. That matters
//   more here than for unordered_dense: absl uses a custom hash functor
//   verbatim, taking the low bits for the slot index and bits 57-63 for the
//   control byte. A weakly-mixed hash would collide in both at once. XXH3 is
//   fine; do not replace it with a raw bit-packing of the k-mer.
// - prefetch() is a no-op, as in the unordered_dense table. absl does expose
//   prefetch_hash(), but it warms the control bytes for a hash we would have to
//   compute early, and count.hpp's prefetch hook is handed a Kmer_Window rather
//   than a key. Left for later if the measurements are close.
// - Single-threaded per partition only, like the other two with mt_=false.

#include "kache-hash/Kmer.hpp"
#include "kache-hash/Streaming_Kmer_Hash_Table.hpp"   // ResizeEvent, Kmer_Window

#include <absl/container/flat_hash_map.h>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace tuna {

// to_u64() is XXH3: already a well-mixed 64-bit hash, handed to absl as-is.
template <uint16_t k>
struct Absl_Kmer_Hash
{
    uint64_t operator()(const kache_hash::Kmer<k>& key) const noexcept { return key.to_u64(); }
};

template <uint16_t k, bool mt_, typename T_ = uint32_t, uint16_t l = 19, bool canonical_ = true>
class Absl_Kmer_Hash_Table
{
    static_assert(!mt_, "Absl_Kmer_Hash_Table is single-threaded; one table per partition.");

public:
    using key_t = kache_hash::Kmer<k>;
    using val_t = T_;
    using map_t = absl::flat_hash_map<key_t, val_t, Absl_Kmer_Hash<k>>;

    // Empty for the same reason as the unordered_dense table: nothing is
    // shared, and the token only keeps the call sites identical.
    struct Token {};

    // Label for an always-empty overflow report; see count.hpp.
    static constexpr std::size_t OV_HIST_BITS_PUBLIC = 8;

private:
    map_t map_;

public:
    // `resize_worker_c` and `lf` are accepted and ignored: absl manages its own
    // growth at a fixed 7/8 max load factor.
    explicit Absl_Kmer_Hash_Table(std::size_t max_sz = 1 << 22,
                                  uint64_t resize_worker_c = 1,
                                  double lf = 0.8)
    {
        (void)resize_worker_c; (void)lf;
        map_.reserve(max_sz);
    }

    Absl_Kmer_Hash_Table(const Absl_Kmer_Hash_Table&) = delete;
    Absl_Kmer_Hash_Table& operator=(const Absl_Kmer_Hash_Table&) = delete;

    // `f` is a transform, not a mutator: an existing count becomes f(count) and
    // a fresh key is stored as `val`. Returns the count before the update, 0 on
    // first insertion, matching the kache table.
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

    void prefetch(const kache_hash::Kmer_Window<k, l>&) const noexcept {}
    void prefetch_packed(const uint8_t*, uint16_t) const noexcept {}

    template <typename F>
    void for_each(F&& f) const { for (const auto& e : map_) f(e); }

    std::size_t size() const { return map_.size(); }
    std::size_t capacity() const { return map_.capacity(); }
    std::size_t bucket_count() const { return map_.bucket_count(); }

    // No overflow table and no resize log; count.hpp reads these
    // unconditionally, so they answer "nothing here".
    uint64_t overflow_insert_count() const { return 0; }
    std::vector<kache_hash::ResizeEvent> resize_log() const { return {}; }
    std::vector<std::pair<uint32_t, uint64_t>> overflow_top_minimizers(std::size_t) const { return {}; }
};

} // namespace tuna
