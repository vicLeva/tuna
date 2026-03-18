
#ifndef KACHE_HASH_STREAMING_KMER_HASHSET_HPP
#define KACHE_HASH_STREAMING_KMER_HASHSET_HPP



#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "Kmer_Utility.hpp"
#include "Directed_Vertex.hpp"
#include "Rolling_Hash.hpp"
#include "Concurrent_Hash_Table.hpp"
#include "RW_Lock.hpp"
#include "kache-hash/DNA.hpp"
#include "utility.hpp"
#include "minimizer_window.hpp"

#include <cstdint>
#include <cstddef>
#include <bit>
#include <atomic>
#include <optional>
#include <utility>
#include <type_traits>
#include <cstring>
#include <iostream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <cassert>

#include <immintrin.h>


namespace kache_hash
{


#define KACHE_HASH_LIKELY(cond)     __builtin_expect(cond, 1)
#define KACHE_HASH_UNLIKELY(cond)   __builtin_expect(cond, 0)


template <uint16_t k, uint16_t l> class Kmer_Window;


// Hash table particularly suited for streaming `k`-mer operations. `mt_`
// denotes whether a concurrent table is desired. If a hash map is required, the
// mapped type is `T_`. `l`-minimizers are used in underlying computations.
// NB: concurrent insertions are guaranteed to function correctly only if
// deletions are not performed, as of now.
// Per-resize diagnostic record populated by Streaming_Kmer_Hash_Table::try_resize().
struct ResizeEvent {
    bool     overflow_triggered; // true = ov_approx_sz_ hit ov_resize_th_; false = main table load
    double   elapsed_s;          // wall-clock duration of resize()
    uint64_t old_cap;            // capacity (k-mer slots) before doubling
    uint64_t ov_count;           // ov_approx_sz_ snapshot at trigger time
    uint64_t main_sz;            // approx_sz snapshot at trigger time
};


template <uint16_t k, bool mt_, typename T_ = void, uint16_t l = 19>
class Streaming_Kmer_Hash_Table
{
private:

    static constexpr double lf_default = 0.8;   // Default load factor: 80%.
    static constexpr double of_default = 0.3;   // Fixed overflow factor: 30% — raises ov_resize_th_ from 5% to 15% of capacity, reducing cascade resizes for hot minimizers.
    static constexpr std::size_t B = 32;    // Bucket size.
    static_assert(B % 8 == 0, "Bucket size must be a multiple of 8.");
    static_assert(B == 32, "Bucket size needs to be 32 for the current vectorization scheme to function.");

    static constexpr bool is_set_ = std::is_void_v<T_>; // Whether a hash-set.
    static constexpr bool is_map_ = !is_set_;   // Whether a hash-map.
    typedef std::conditional_t<is_map_, std::pair<Kmer<k>, T_>, Kmer<k>> flat_t;    // Underlying data type in the flat table.
    typedef std::optional<T_> val_t;    // Type of values returned by hash-map operations.
    typedef std::conditional_t<is_map_, val_t, bool> find_ret_t;    // Return type for find / insert methods.
    static constexpr auto null_val = [](){ if constexpr(is_map_) return std::nullopt;   else return false; }(); // Marker for absent keys (also values).

    static constexpr auto& key(const flat_t& e) { if constexpr(is_map_) return e.first; else return e; }    // Returns the key from a table-element `e`.
    static constexpr auto& val(const flat_t& e) { if constexpr(is_map_) return e.second; else return e; }   // Returns the value from a table-element `e`.

    static constexpr uint8_t min_orientation_mask = 0b0100'0000;

    // Seed used to derive bucket routing hash from the minimizer's ntHash value.
    // Distinct from the cuckoo probe seeds (0 and 1<<63) used for secondary/tertiary buckets.
    // ntHash is used only for minimizer *selection* (min comparison); bucket placement
    // uses XXH3(ntHash, kBucketSeed) so that the two concerns are decoupled.
    static constexpr uint64_t kBucketSeed = 0x9e3779b97f4a7c15ULL;

    const double lf;    // Load factor.

    std::size_t cap_;   // Capacity of the table in terms of buckets; adjusted to be a power of 2.
    uint64_t resize_th; // At what size of the table to trigger the next resize.

    // Metadata of the k-mers in a bucket.
    struct Metadata
    {
        uint8_t cs[B];  // Checksums of the k-mers.
        uint8_t min_coord[B];   // Minimizer coordinates within the k-mers.

        static constexpr uint8_t lock_mask = 0b1000'0000;   // Bit of `min_coord[0]` used to lock the bucket.
    };

    struct Overflow_Set_Val
    {
        uint8_t cs;
        uint8_t min_coord;

        Overflow_Set_Val() {}
        Overflow_Set_Val(uint8_t cs, uint8_t min_coord): cs(cs), min_coord(min_coord) {}
    };

    struct Overflow_Map_Val
    {
        uint8_t cs;
        uint8_t min_coord;
        T_ val;

        Overflow_Map_Val() {}
        Overflow_Map_Val(uint8_t cs, uint8_t min_coord, const T_& val): cs(cs), min_coord(min_coord), val(val) {}
    };

    typedef std::conditional_t<is_set_,
        Concurrent_Hash_Table<Kmer<k>, Overflow_Set_Val, Hash<Kmer<k>>>,
        Concurrent_Hash_Table<Kmer<k>, Overflow_Map_Val, Hash<Kmer<k>>> > ov_t; // The overflow table type.


    flat_t* T;  // Flat table of size `cap x B`.
    Metadata* M;    // Metadata table of the buckets.
    ov_t* ov;   // The overflow table.

    flat_t* T_new;  // New flat table after a resize.
    Metadata* M_new;    // New metadata table of buckets after a resize.
    ov_t* ov_new;    // New overflow table after a resize.

    std::atomic_uint16_t registered_thread_count;   // Number of threads registered to use the table.
    std::atomic_uint64_t approx_sz; // Approximate size (lower bound) of the hash table.
    std::atomic_uint64_t ov_approx_sz_;    // Overflow inserts since last resize (drives ov_resize_th_; resets on resize).
    std::atomic_uint64_t ov_total_inserts_; // Cumulative overflow inserts across all resizes (never reset).
    uint64_t ov_resize_th_;                // Overflow insert count at which to trigger resize.

    // Histogram of overflow inserts by minimizer hash (top 16 bits of nt_h → 65536 bins).
    // Populated by the public upsert when a new k-mer goes to overflow.
    static constexpr std::size_t OV_HIST_BITS = 8;  // 256 bins × 8 B = 2 KB — fully L1-resident, eliminating LLC misses per overflow insert (vs 512 KB at 16 bits).
    static constexpr std::size_t OV_HIST_SIZE = std::size_t(1) << OV_HIST_BITS;
    std::vector<std::atomic_uint64_t> ov_min_hist_;

    // Thread-local flag: set by the overflow insert path, read by the public upsert.
    static thread_local bool tl_ov_happened_;
    std::vector<Padded<uint64_t>> checkpoint;   // `checkpoint[i]` is the number of elements remaining to be added by the `i`'th user thread before it checks for a resize.
    static constexpr uint64_t resize_checkpoint = 16384;    // Checkpoint number of elements for a user thread to add before checking for resize.

    static constexpr uint64_t max_user = 64;    // Maximum number of users allowed to use the table.
    RW_Lock<max_user> table_lock; // Readers-writers lock to resize the table exclusively.

    const uint64_t resize_worker_c; // Number of worker threads to spawn for a concurrent resize.

    std::vector<ResizeEvent> resize_log_; // Per-table resize event log (populated by try_resize).


    // Returns `true` iff `key` matches to the key of the element `e`.
    static constexpr bool key_equals(const flat_t& e, const Kmer<k>& key) { return Streaming_Kmer_Hash_Table::key(e) == key; }

    // Returns the 32-bit equality mask of the vectors `x` and `y`.
    static uint32_t eq_mask(__m256i x, __m256i y);

    // Returns the size of the `b`'th bucket.
    std::size_t bucket_size(std::size_t b) const;

    // Returns the size of the main table.
    std::size_t main_table_size() const;

    // Locks the `b`'th bucket for exclusive-write access in concurrent
    // settings. `resize_` denotes whether the bucket requires to be locked for
    // resizing or not.
    template <bool resize_ = false> void lock(std::size_t b);

    // Locks the `p`'th and the `q`'th buckets for exclusive-write accesses in
    // concurrent settings, taking care of circular deadlocks and when the
    // buckets are the same. `resize_` denotes whether the buckets require to be
    // locked for resizing or not.
    template <bool resize_ = false> void lock_ordered(std::size_t p, std::size_t q);

    // Unlocks the `b`'th bucket from exclusive-write access in concurrent
    // settings. `resize_` denotes whether the bucket requires to be unlocked
    // during resizing or not.
    template <bool resize_ = false> void unlock(std::size_t b);

    // Unlocks the `p` and the `q`'th buckets from exclusive-write access in
    // concurrent settings, taking care of circular deadlocks and when the
    // buckets are the same. `resize_` denotes whether the buckets require to be
    // unlocked during resizing or not.
    template <bool resize_ = false> void unlock_ordered(std::size_t p, std::size_t q);

    // Returns the associated value iff the k-mer `key` with checksum `c` is
    // present in the `i`'th bucket.
    find_ret_t find(Kmer<k> key, uint8_t c, std::size_t i) const;

    // Returns the associated value iff the k-mer `key` is present in the
    // overflow table.
    find_ret_t find_in_overflow(Kmer<k> key) const;

    // Updates the associated value of the k-mer `key` with checksum `c` in the
    // `i`'th bucket through `f` iff it is present there, returning the earlier
    // value. Returns `null_opt` otherwise.
    template <typename F, typename V = T_>
    val_t find_and_update(Kmer<k> key, uint8_t c, std::size_t i, const F& f) requires (is_map_);

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j`.
    void place_at(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, const std::size_t j)
    { T[i * B + j] = e, M[i].cs[j] = c, M[i].min_coord[j] = m; }

    // Directly places the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket at index `j` of the new resized
    // table.
    void place_at_resize(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, const std::size_t j)
    { T_new[i * B + j] = e, M_new[i].cs[j] = c, M_new[i].min_coord[j] = m; }

    // Tries to insert the element `e` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket. It succeeds iff the key was absent
    // in the bucket beforehand and there was empty space, returning `true` in
    // that case. Returns `false` otherwise; in this case if the key was already
    // in the bucket, it puts the associated value at `val`, or else puts
    // `null_val` at `val`.
    bool try_insert_at(const flat_t& e, uint8_t c, uint8_t m, std::size_t i, find_ret_t& val);

    // Tries to insert the element `x` with checksum `c` and minimizer-
    // coordinate `m` into the `i`'th bucket of the newly resized table. It
    // succeeds iff there was empty space, returning `true` in that case.
    // Returns `false` otherwise.
    bool try_insert_at_resize(const flat_t& x, uint8_t c, uint8_t m, std::size_t i);

    // Tries to update the associated value of the k-mer `key` with checksum `c`
    // and minimizer-coordinate `m` at the `i`'th bucket through `f` iff it is
    // present there, and also puts the associated value to `old_val`. Otherwise
    // if the bucket has space, puts `key` with value `val` into it and puts
    // `null_val` at `old_val`. Returns `true` iff the upsert succeeds.
    template <typename F, typename V = T_>
    bool try_upsert_at(Kmer<k> key, uint8_t c, uint8_t m, std::size_t i, const F& f, const V& val, val_t& old_val) requires (is_map_);

    // Tries to insert the element `e` with checksum `c`, minimizer-hash `h`,
    // and -coordinate `m` into the table. It succeeds iff the key was absent in
    // the table beforehand, and returns `null_val` in that case. Otherwise it
    // returns the currently associated value.
    find_ret_t insert(const flat_t& e, uint8_t c, uint64_t h, uint8_t m);

    // Inserts the element `x` with checksum `c` and minimizer-coordinate `m`
    // into the resized table.
    void insert_at_resize(const flat_t& x, uint8_t c, uint8_t m);

    // Inserts the k-mer `key` with checksum `c`, minimizer-hash `h`, and
    // -coordinate `m` into the table with associated value `val` iff it is not
    // already present. Otherwise updates the currently associated value through
    // `f`. Returns the old associated value if the key was already present, and
    // returns `null_val` otherwise.
    template <typename F, typename V = T_>
    val_t upsert(const Kmer<k> key, uint8_t c, uint64_t h, uint8_t m, const F& f, const V& val) requires (is_map_);

    // Tries to resize the hash table to double its capacity.
    void try_resize();

    // Resizes the hash table to double its capacity.
    void resize();

public:

    class Token;

    // Constructs a streaming k-mer hash table to support up-to `max_sz` k-mers,
    // with a maximum load-factor of `lf`.
    Streaming_Kmer_Hash_Table(std::size_t max_sz = 1 << 22, uint64_t resize_worker_c = 16, double lf = lf_default);

    // Deserializes a streaming k-mer hash table from the stream `is`.
    // Streaming_Kmer_Hash_Table(std::ifstream& is);

    Streaming_Kmer_Hash_Table(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table(Streaming_Kmer_Hash_Table&&) = delete;
    Streaming_Kmer_Hash_Table& operator=(const Streaming_Kmer_Hash_Table&) = delete;
    Streaming_Kmer_Hash_Table& operator=(Streaming_Kmer_Hash_Table&&) = delete;

    ~Streaming_Kmer_Hash_Table();

    // Number of bits used for the overflow histogram (public so callers can label output).
    static constexpr std::size_t OV_HIST_BITS_PUBLIC = OV_HIST_BITS;

    // Returns the capacity (in terms of keys) of the table.
    auto capacity() const { return cap_ * B; }

    // Returns the number of buckets in the flat table.
    std::size_t bucket_count() const { return cap_; }

    // Returns the size of the table.
    std::size_t size() const;

    // Returns the size of the overflow table.
    std::size_t overflow_size() const { return ov->size();  }

    // Returns the total number of new k-mers inserted into overflow (cumulative across all resizes).
    uint64_t overflow_insert_count() const { return ov_total_inserts_.load(std::memory_order_relaxed); }

    // Returns the top-N (minimizer_bin, overflow_count) pairs sorted descending.
    // minimizer_bin = top OV_HIST_BITS bits of the canonical ntHash minimizer value.
    std::vector<std::pair<uint32_t, uint64_t>> overflow_top_minimizers(std::size_t n) const
    {
        std::vector<std::pair<uint32_t, uint64_t>> result;
        result.reserve(OV_HIST_SIZE);
        for(std::size_t i = 0; i < OV_HIST_SIZE; ++i)
        {
            const auto cnt = ov_min_hist_[i].load(std::memory_order_relaxed);
            if(cnt > 0) result.emplace_back(static_cast<uint32_t>(i), cnt);
        }
        std::partial_sort(result.begin(),
                          result.begin() + std::min(n, result.size()),
                          result.end(),
                          [](const auto& a, const auto& b){ return a.second > b.second; });
        if(result.size() > n) result.resize(n);
        return result;
    }

    // Returns the log of resize events recorded during insertion.
    const std::vector<ResizeEvent>& resize_log() const { return resize_log_; }

    // Clears the table.
    void clear();

    // Inserts the k-mer defined by the window `w` into the table iff it is not
    // already present in the table, returning `false` in that case. Returns
    // `true` otherwise; i.e. returns the presence-status before the insertion.
    // The invoking user thread must provide its token `token`.
    bool insert(const Kmer_Window<k, l>& w, const Token& token) requires (is_set_);

    // Inserts the k-mer defined by the window `w` and associated value `val`
    // into the table iff the key is not already present in the table, returning
    // `null_opt` in that case. Returns the associated current value otherwise;
    // i.e. returns the associated value pre-insert. The invoking user thread
    // must provide its token `token`.
    template <typename V = T_> val_t insert(const Kmer_Window<k, l>& w, const V& val, const Token& token) requires (is_map_);

    // Updates the associated value of the k-mer defined by `w` in the table
    // through `f`, iff it is present there. Otherwise inserts it into the table
    // with associated value `val`. Returns the old associated value if the key
    // was already present in the table, and returns `null_val` otherwise. The
    // invoking user thread must provide its token `token`.
    template <typename F, typename V = T_>
    val_t upsert(const Kmer_Window<k, l>& w, const F& f, const V& val, const Token& token) requires (is_map_);

    // Returns `true` iff the k-mer defined by the window `w` is present in the
    // table.
    find_ret_t find(const Kmer_Window<k, l>& w) const;

    // Issues a non-blocking L1 prefetch for the primary bucket of the k-mer
    // defined by `w`.  Call this PREFETCH_DIST iterations before the matching
    // upsert to hide LLC miss latency (~200 ns) behind useful computation.
    void prefetch(const Kmer_Window<k, l>& w) const
    {
        // Prefetch the metadata cache line (cs[32] + min_coord[32] = 64 bytes).
        // The SIMD checksum scan always hits M[b] first; T[b*B+j] is only
        // touched on a checksum match, which is rare at low load factors.
        const auto pf_nt_h = w.minimizer_hash();
        const auto pf_h = XXH3_64bits_withSeed(&pf_nt_h, sizeof(pf_nt_h), kBucketSeed);
        __builtin_prefetch(&M[pf_h & (cap_ - 1)], 0, 3);
    }

    // Issues a non-blocking L1 prefetch using only the l-mer minimizer stored at
    // offset `min_pos` inside a packed superkmer byte array (kache 2-bit encoding).
    // O(l) work — cheaper than prefetch(Kmer_Window) which needs a fully
    // initialised window.  Call this one superkmer ahead to hide the LLC miss
    // behind the current superkmer's hot loop.
    void prefetch_packed(const uint8_t* packed, const uint8_t min_pos) const
    {
        static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
        char buf_l[l];
        for (uint16_t i = 0; i < l; ++i) {
            const uint16_t pos = static_cast<uint16_t>(min_pos) + i;
            buf_l[i] = B2C[(packed[pos >> 2] >> (6u - 2u * (pos & 3u))) & 3u];
        }
        nt_hash::Roller<l> roller;
        roller.init(buf_l);
        const uint64_t pf_nt_h = roller.canonical();
        const uint64_t pf_h = XXH3_64bits_withSeed(&pf_nt_h, sizeof(pf_nt_h), kBucketSeed);
        __builtin_prefetch(&M[pf_h & (cap_ - 1)], 0, 3);
    }

    // Registers a new user and returns a unique token for it.
    Token register_user() { return Token(registered_thread_count++); }

    class Iterator;

    // Returns an iterator to the beginning of the table.
    Iterator cbegin() const { return Iterator(*this, 0, 0); }

    // Returns an iterator to the (exclusive-)end of the table.
    Iterator cend() const { return Iterator(*this, cap_, -1); };

    // Invokes `f` with a `const flat_t&` for every element in the table,
    // including elements residing in the overflow table. For hash-maps,
    // `flat_t` is `std::pair<Kmer<k>, T_>`; for hash-sets it is `Kmer<k>`.
    template <typename F> void for_each(F&& f) const;

    // Serializes the table to the stream `os`.
    // void serialize(std::ofstream& os) const;

    // Deserializes the table from the stream `is`.
    // void deserialize(std::ifstream& is);

    // Reports developer statistics about the table.
    void report_stats() const;
};


// Token required to use the hash table.
template <uint16_t k, bool mt_, typename T_, uint16_t l>
class alignas(L1_CACHE_LINE_SIZE)
    Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Token
{
    friend class Streaming_Kmer_Hash_Table;

private:

    uint16_t id;    // Unique ID of the token.

    Token(const uint16_t id): id(id)
    {}

public:

    Token(): id(0)
    {}
};


// A window defining a k-mer with its associated metadata: its hash and
// minimizer.  The minimizer uses canonical ntHash (MinimizerWindow<k>),
// replacing the earlier XXH3-based Min_Iterator<k,l>.
template <uint16_t k, uint16_t l>
class Kmer_Window
{
    template <uint16_t k_, bool mt_, typename T_, uint16_t l_> friend class Streaming_Kmer_Hash_Table;

    Directed_Vertex<k>    v;
    Rolling_Hash<k, true> rh;
    MinimizerWindow<k, l> nt_min;
    uint64_t              precomp_nt_h_ = 0;
    bool                  use_precomp_  = false;

public:

    // Initializes the k-mer window at the beginning of the ASCII sequence `s`.
    void init(const char* const s)
    {
        v = Directed_Vertex<k>(Kmer<k>(s));
        rh.init(s);
        nt_min.reset(s);
        use_precomp_ = false;
    }

    // Initializes the k-mer window from 2-bit packed DNA in kache encoding
    // (A=0, C=1, G=2, T=3), 4 bases/byte MSB-first.
    // `packed` must contain at least k bases.
    void init_packed(const uint8_t* packed)
    {
        static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
        char buf[k];
        for (uint16_t i = 0; i < k; ++i)
            buf[i] = B2C[(packed[i >> 2] >> (6u - 2u * (i & 3u))) & 3u];
        init(buf);
    }

    // Initializes the k-mer window from ASCII sequence `s` with a precomputed
    // canonical ntHash of the minimizer l-mer (e.g. from MinimizerWindow::hash()).
    // Skips MinimizerWindow::reset() — nt_h is reused for all k-mers in the
    // superkmer since they all share the same minimizer.
    void init_with_hash_ascii(const char* const s, uint64_t nt_h)
    {
        v = Directed_Vertex<k>(Kmer<k>(s));
        rh.init(s);
        precomp_nt_h_ = nt_h;
        use_precomp_  = true;
    }

    // Initializes the k-mer window from packed data with a precomputed canonical
    // ntHash of the minimizer l-mer (from the stored min_pos header byte).
    // Skips MinimizerWindow::reset() — nt_h is reused for all k-mers in the
    // superkmer since they all share the same minimizer.
    void init_packed_with_hash(const uint8_t* packed, uint64_t nt_h)
    {
        static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
        char buf[k];
        for (uint16_t i = 0; i < k; ++i)
            buf[i] = B2C[(packed[i >> 2] >> (6u - 2u * (i & 3u))) & 3u];
        v = Directed_Vertex<k>(Kmer<k>(buf));
        rh.init(buf);
        precomp_nt_h_ = nt_h;
        use_precomp_  = true;
    }

    // Fused init: decode k packed bases once, derive ntHash of the l-mer at
    // min_pos from the decoded buffer, init Directed_Vertex + Rolling_Hash.
    // Returns the canonical ntHash of the minimizer l-mer.
    // Replaces the two-call pattern: lmer_nt_hash(packed, min_pos) then
    // init_packed_with_hash(packed, mh) — eliminates one full decode pass.
    uint64_t init_packed_with_min(const uint8_t* packed, uint8_t min_pos)
    {
        static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
        char buf[k];
        for (uint16_t i = 0; i < k; ++i)
            buf[i] = B2C[(packed[i >> 2] >> (6u - 2u * (i & 3u))) & 3u];
        v = Directed_Vertex<k>(Kmer<k>(buf));
        rh.init(buf);
        nt_hash::Roller<l> roller;
        roller.init(buf + min_pos);
        precomp_nt_h_ = roller.canonical();
        use_precomp_  = true;
        return precomp_nt_h_;
    }

    // Initializes the k-mer window from a packed super-kmer word array.
    void init(const uint64_t* super_kmer, const std::size_t word_count)
    {
        v.from_super_kmer(super_kmer, word_count);
        rh.init(v.kmer());
        // Decode kmer to ASCII for nt_min (front-to-back order).
        static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
        char buf[k];
        for (uint16_t i = 0; i < k; ++i)
            buf[i] = B2C[v.kmer().base_at(k - 1 - i)];
        nt_min.reset(buf);
        use_precomp_ = false;
    }

    // Advances the window by one nucleobase `b` (kache encoding: A=0,C=1,G=2,T=3).
    void advance(const DNA::Base b)
    {
        v.roll_forward(b);
        if (!use_precomp_) nt_min.advance_kache(static_cast<uint8_t>(b));
        rh.advance(b);
    }

    // Advances the window by one ASCII character `ch`.
    void advance(const char ch)
    {
        assert(!DNA_Utility::is_placeholder(ch));
        const DNA::Base b = DNA_Utility::map_base(ch);
        v.roll_forward(b);
        if (!use_precomp_) nt_min.advance_kache(static_cast<uint8_t>(b));
        rh.advance(b);
    }

    // Returns the ntHash of the l-minimizer of the current k-mer.
    // Returns precomputed hash if provided at init.
    auto minimizer_hash() const { return use_precomp_ ? precomp_nt_h_ : nt_min.hash(); }

    // Same as minimizer_hash() with min_coord output (used by upsert/insert).
    uint64_t minimizer_nt_hash(uint8_t& m) const
    {
        if (use_precomp_) { m = 0; return precomp_nt_h_; }
        return nt_min.hash(m);
    }

    // Returns the hash of the current k-mer in the forward-strand.
    auto hash_fwd() const { return rh.hash_fwd(); }

    // Returns the hash of the current k-mer in the reverse-strand.
    auto hash_rev() const { return rh.hash_rev(); }
};


template <uint16_t k, bool mt_, typename T_, uint16_t l>
class Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Iterator
{
    friend class Streaming_Kmer_Hash_Table<k, mt_, T_, l>;
    typedef Streaming_Kmer_Hash_Table<k, mt_, T_, l> ht_t;

private:

    const ht_t& ht; // Associated hash table.
    std::size_t b;  // Index of the bucket in the table where the iterator sits.
    int8_t slot;    // Index of the slot in the bucket where the iterator sits.

    // Constructs an iterator for the hash table `ht` that sits in bucket `b` at
    // slot `s`.
    Iterator(const ht_t& ht, const std::size_t b, const int8_t s): ht(ht), b(b), slot(s)
    {
        assert(b == 0 || b == ht.cap_);
        if(b < ht.cap_ && s >= 0 && s < static_cast<int8_t>(B) && !ht.M[b].cs[s])
            operator++();
    }

public:

    // Repositions the iterator to the next successive element in the map, or to
    // its end if there are no more elements.
    void operator++()
    {
        for(; b < ht.cap_; ++b)
        {
            for(slot++; slot < static_cast<int8_t>(B); ++slot)
                if(ht.M[b].cs[slot])
                    return;

            slot = -1;
        }
    }

    const ht_t::flat_t& operator*() { assert(b < ht.cap_ && 0 <= slot && slot < static_cast<int8_t>(B)); return ht.T[b * B + slot]; }

    auto checksum() const { assert(b < ht.cap_ && 0 <= slot && slot < static_cast<int8_t>(B)); return ht.M[b].cs[slot]; }

    auto min_coord() const { assert(b < ht.cap_ && 0 <= slot && slot < static_cast<int8_t>(B)); return ht.M[b].min_coord[slot]; }

    bool operator==(const Iterator& rhs) const { return &ht == &rhs.ht && b == rhs.b && slot == rhs.slot; }
};


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::Streaming_Kmer_Hash_Table(const std::size_t max_sz, const uint64_t resize_worker_c, const double lf):
      lf(lf)
    , cap_(ceil_pow_2((max_sz / lf) / B))
    , resize_th(capacity() * lf)
    , T(aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE))
    , M(aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE))
    , ov(new ov_t(capacity() * of_default))
    , T_new(nullptr)
    , M_new(nullptr)
    , ov_new(nullptr)
    , registered_thread_count(0)
    , approx_sz(0)
    , ov_approx_sz_(0)
    , ov_total_inserts_(0)
    , ov_resize_th_(static_cast<uint64_t>(capacity() * of_default * 0.5))
    , ov_min_hist_(OV_HIST_SIZE)
    , checkpoint(std::thread::hardware_concurrency(), resize_checkpoint)
    , resize_worker_c(resize_worker_c)
{
    clear();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline Streaming_Kmer_Hash_Table<k, mt_, T_, l>::~Streaming_Kmer_Hash_Table()
{
    kache_hash::deallocate(T), kache_hash::deallocate(M), delete(ov);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::clear()
{
    std::memset(M, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    approx_sz = 0;
    ov_approx_sz_.store(0, std::memory_order_relaxed);
    for(auto& bin : ov_min_hist_) bin.store(0, std::memory_order_relaxed);
    std::for_each(checkpoint.begin(), checkpoint.end(), [](auto& c){ c = resize_checkpoint; });

    ov->clear();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline uint32_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::eq_mask(const __m256i x, const __m256i y)
{
#ifdef AVX_512
    return _mm256_cmpeq_epi8_mask(x, y);
#else
    const auto cmp = _mm256_cmpeq_epi8(x, y);
    return static_cast<uint32_t>(_mm256_movemask_epi8(cmp));
#endif
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::bucket_size(const std::size_t b) const
{
    const auto zero  = _mm256_setzero_si256();
    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[b].cs));
    const auto cmp   = eq_mask(c_vec, zero);

    assert(std::countr_zero(cmp) + std::countl_one(cmp) == B);
    return std::countr_zero(cmp);   // This must hold if no deletion is performed at the table.
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::main_table_size() const
{
    std::size_t sz = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        sz += bucket_size(b);

    return sz;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock(const std::size_t b)
{
    if constexpr(mt_)
        while(true)
        {
            if constexpr(!resize_)
            {
                const auto exp_unlocked = M[b].min_coord[0] & ~Metadata::lock_mask;
                if(__sync_bool_compare_and_swap(&M[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
                    break;
            }
            else
            {
                const auto exp_unlocked = M_new[b].min_coord[0] & ~Metadata::lock_mask;
                if(__sync_bool_compare_and_swap(&M_new[b].min_coord[0], exp_unlocked, exp_unlocked | Metadata::lock_mask))
                    break;
            }
        }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::lock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_)
        KACHE_HASH_LIKELY(p != q) ?
            lock<resize_>(std::min(p, q)), lock<resize_>(std::max(p, q)) :
            lock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock(const std::size_t b)
{
    if constexpr(mt_)
    {
        if constexpr(!resize_)
        {
            assert(M[b].min_coord[0] & Metadata::lock_mask);
            M[b].min_coord[0] &= ~Metadata::lock_mask;
        }
        else
        {
            assert(M_new[b].min_coord[0] & Metadata::lock_mask);
            M_new[b].min_coord[0] &= ~Metadata::lock_mask;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <bool resize_>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::unlock_ordered(const std::size_t p, const std::size_t q)
{
    if constexpr(mt_)
        KACHE_HASH_LIKELY(p != q) ?
                unlock<resize_>(std::max(p, q)), unlock<resize_>(std::min(p, q)) :
                unlock<resize_>(p);
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline std::size_t Streaming_Kmer_Hash_Table<k, mt_, T_, l>::size() const
{
    // return main_table_size() + overflow_size();
    std::size_t sz = approx_sz;
    std::for_each(checkpoint.cbegin(), checkpoint.cend(), [&](auto& c){ sz += c.unwrap() > 0 ? resize_checkpoint - c.unwrap() : 0; });
    return sz;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const Kmer<k> key, const uint8_t c, const std::size_t i) const -> find_ret_t
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    const auto cmp = eq_mask(c_vec, query);

    const auto j = cmp ? __builtin_ctz(cmp) : 0;
    if(M[i].cs[j] && key_equals(T[i * B + j], key))
    {
        if constexpr(is_map_)   return T[i * B + j].second;
        else                    return true;
    }

    auto cmp_rem = cmp & (cmp - 1);
    while(cmp_rem)
    {
        const auto j = __builtin_ctz(cmp_rem);
        if(key_equals(T[i * B + j], key))
        {
            if constexpr(is_map_)   return T[i * B + j].second;
            else                    return true;
        }

        cmp_rem &= (cmp_rem - 1);
    }

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find_in_overflow(const Kmer<k> key) const -> find_ret_t
{
    if constexpr(is_set_)
        return ov->find(key);
    else
    {
        const auto p = ov->find(key);
        return !p ? null_val : val_t(p->val);
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find_and_update(const Kmer<k> key, const uint8_t c, const std::size_t i, const F& f) -> val_t requires (is_map_)
{
    assert(i < cap_);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto query = _mm256_set1_epi8(static_cast<char>(c));
    auto cmp = eq_mask(c_vec, query);

    while(cmp)
    {
        const auto j = __builtin_ctz(cmp);
        if(key_equals(T[i * B + j], key))
        {
            const auto prev_val = T[i * B + j].second;
            T[i * B + j].second = f(prev_val);
            return prev_val;
        }

        cmp &= (cmp - 1);
    }

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i, find_ret_t& val)
{
    lock(i);

    const auto& key = this->key(e);
    if(const auto r = find(key, c, i))
    {
        unlock(i);
        val = r;
        return false;
    }

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock(i);
        val = null_val;
        return false;
    }

    place_at(e, c, m, i, j);
    unlock(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_insert_at_resize(const flat_t& e, const uint8_t c, const uint8_t m, const std::size_t i)
{
    lock<true>(i);

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M_new[i].cs));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock<true>(i);
        return false;
    }

    place_at_resize(e, c, m, i, j);
    unlock<true>(i);
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_upsert_at(const Kmer<k> key, const uint8_t c, const uint8_t m, const std::size_t i, const F& f, const V& val, val_t& old_val) requires (is_map_)
{
    lock(i);

    if(const auto r = find_and_update(key, c, i, f))
    {
        unlock(i);
        old_val = r;
        return true;
    }

    const auto c_vec = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[i].cs));
    const auto zero = _mm256_setzero_si256();
    const auto cmp = eq_mask(c_vec, zero);

    const auto j = std::countr_zero(cmp);
    if(j == B)
    {
        unlock(i);
        return false;
    }

    place_at(std::make_pair(key, val), c, m, i, j);
    unlock(i);

    old_val = null_val;
    return true;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const Kmer_Window<k, l>& w, const Token& token) requires (is_set_)
{
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto nt_h = w.minimizer_nt_hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto h = XXH3_64bits_withSeed(&nt_h, sizeof(nt_h), kBucketSeed);

    if constexpr(mt_)   table_lock.lock_shared(token.id);
    const auto r = insert(w.v.canonical(), c, h, m);
    if constexpr(mt_)   table_lock.unlock_shared(token.id);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        if(approx_sz >= resize_th || ov_approx_sz_.load(std::memory_order_relaxed) >= ov_resize_th_)
            try_resize();
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const Kmer_Window<k, l>& w, const V& val, const Token& token) -> val_t requires (is_map_)
{
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto nt_h = w.minimizer_nt_hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto h = XXH3_64bits_withSeed(&nt_h, sizeof(nt_h), kBucketSeed);

    if constexpr(mt_)   table_lock.lock_shared(token.id);
    const auto r = insert(std::make_pair(w.v.canonical(), val), c, h, m);
    if constexpr(mt_)   table_lock.unlock_shared(token.id);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        if(approx_sz >= resize_th || ov_approx_sz_.load(std::memory_order_relaxed) >= ov_resize_th_)
            try_resize();
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::upsert(const Kmer_Window<k, l>& w, const F& f, const V& val, const Token& token) -> val_t requires (is_map_)
{
    (void)token;
    uint8_t m;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));  // Avoiding checksum 0 by overloading checksum 1.
    const auto nt_h = w.minimizer_nt_hash(m);
    m = w.v.in_canonical_form() ? m : (m ^ min_orientation_mask);
    const auto h = XXH3_64bits_withSeed(&nt_h, sizeof(nt_h), kBucketSeed);

    tl_ov_happened_ = false;
    if constexpr(mt_)   table_lock.lock_shared(token.id);
    const auto r = upsert(w.v.canonical(), c, h, m, f, val);
    if constexpr(mt_)   table_lock.unlock_shared(token.id);

    if(tl_ov_happened_)
        ov_min_hist_[nt_h >> (64 - OV_HIST_BITS)].fetch_add(1, std::memory_order_relaxed);

    auto& cp = checkpoint[token.id].unwrap();
    cp -= (r == null_val);
    if(KACHE_HASH_UNLIKELY(cp == 0))
    {
        approx_sz += resize_checkpoint;
        cp = resize_checkpoint;
        if(approx_sz >= resize_th || ov_approx_sz_.load(std::memory_order_relaxed) >= ov_resize_th_)
            try_resize();
    }

    return r;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert(const flat_t& x, const uint8_t c, const uint64_t h, const uint8_t m) -> find_ret_t
{
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    find_ret_t v;
    if(const auto r = try_insert_at(x, c, m, b, v))
        return null_val;
    else if(v)
        return v;
    else
    {
        // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
        const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
        const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

        const auto s_1 = h_1 & idx_mask;
        const auto s_2 = h_2 & idx_mask;
        const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);

        const auto key = this->key(x);

        lock_ordered(p, q);
        if(const auto r_p = find(key, c, p))
        {
            unlock_ordered(p, q);
            return r_p;
        }
        else if(const auto r_q = find(key, c, q))
        {
            unlock_ordered(p, q);
            return r_q;
        }
        else
        {
            const auto zero = _mm256_setzero_si256();

            const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
            const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

            const auto cmp_p = eq_mask(c_vec_p, zero);
            const auto cmp_q = eq_mask(c_vec_q, zero);

            const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
            const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
            assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

            if(sz_p < sz_q)
            {
                place_at(x, c, m, p, sz_p);
                unlock_ordered(p, q);
            }
            else if(sz_q < B)
            {
                place_at(x, c, m, q, sz_q);
                unlock_ordered(p, q);
            }
            else
            {
                unlock_ordered(p, q);
                if constexpr(is_set_)
                {
                    const auto ov_r = ov->insert(x, Overflow_Set_Val(c, m));
                    ov_approx_sz_.fetch_add(!ov_r, std::memory_order_relaxed);
                    ov_total_inserts_.fetch_add(!ov_r, std::memory_order_relaxed);
                    return ov_r;
                }
                else
                {
                    const auto p = ov->insert(x.first, Overflow_Map_Val(c, m, x.second));
                    ov_approx_sz_.fetch_add(!p, std::memory_order_relaxed);
                    ov_total_inserts_.fetch_add(!p, std::memory_order_relaxed);
                    return !p ? null_val : val_t(p->val);
                }
            }

            return null_val;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::insert_at_resize(const flat_t& x, const uint8_t c, const uint8_t m)
{
    const auto idx_mask = cap_ - 1;

    const auto key = this->key(x);
    // Recompute minimizer hash via MinimizerWindow<k> (canonical ntHash).
    static constexpr char B2C[4] = {'A', 'C', 'G', 'T'};
    char buf[k];
    for (uint16_t i = 0; i < k; ++i)
        buf[i] = B2C[key.base_at(k - 1 - i)];
    static thread_local MinimizerWindow<k, l> tmp_win;
    tmp_win.reset(buf);
    const auto nt_h = tmp_win.hash();
    const auto h = XXH3_64bits_withSeed(&nt_h, sizeof(nt_h), kBucketSeed);

    const auto b = h & idx_mask;
    if(try_insert_at_resize(x, c, m, b))
        return;

    const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
    const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);
    const auto s_1 = h_1 & idx_mask;
    const auto s_2 = h_2 & idx_mask;
    const auto p = std::min(s_1, s_2), q = std::max(s_1, s_2);

    lock_ordered<true>(p, q);

    const auto zero = _mm256_setzero_si256();

    const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M_new[p].cs));
    const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M_new[q].cs));

    const auto cmp_p = eq_mask(c_vec_p, zero);
    const auto cmp_q = eq_mask(c_vec_q, zero);

    const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
    const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
    assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

    if(sz_p < sz_q)
    {
        place_at_resize(x, c, m, p, sz_p);
        unlock_ordered<true>(p, q);
    }
    else if(sz_q < B)
    {
        place_at_resize(x, c, m, q, sz_q);
        unlock_ordered<true>(p, q);
    }
    else
    {
        unlock_ordered<true>(p, q);
        if constexpr(is_set_)
            ov_new->insert(x, Overflow_Set_Val(c, m));
        else
            ov_new->insert(x.first, Overflow_Map_Val(c, m, x.second));
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F, typename V>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::upsert(const Kmer<k> key, const uint8_t c, const uint64_t h, const uint8_t m, const F& f, const V& val) -> val_t requires (is_map_)
{
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    val_t old_val;
    if(try_upsert_at(key, c, m, b, f, val, old_val))
        return old_val;
    else
    {
        // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
        const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
        const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

        const auto p = h_1 & idx_mask;
        const auto q = h_2 & idx_mask;

        lock_ordered(p, q);
        if(const auto r_p = find_and_update(key, c, p, f))
        {
            unlock_ordered(p, q);
            return r_p;
        }
        else if(const auto r_q = find_and_update(key, c, q, f))
        {
            unlock_ordered(p, q);
            return r_q;
        }
        else
        {
            const auto zero = _mm256_setzero_si256();

            const auto c_vec_p = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[p].cs));
            const auto c_vec_q = _mm256_load_si256(reinterpret_cast<const __m256i*>(M[q].cs));

            const auto cmp_p = eq_mask(c_vec_p, zero);
            const auto cmp_q = eq_mask(c_vec_q, zero);

            const auto sz_p = static_cast<std::size_t>(std::countr_zero(cmp_p));
            const auto sz_q = static_cast<std::size_t>(std::countr_zero(cmp_q));
            assert(std::countr_zero(cmp_p) + std::countl_one(cmp_p) == B), assert(std::countr_zero(cmp_q) + std::countl_one(cmp_q) == B);

            if(sz_p < sz_q)
            {
                place_at(std::make_pair(key, val), c, m, p, sz_p);
                unlock_ordered(p, q);
            }
            else if(sz_q < B)
            {
                place_at(std::make_pair(key, val), c, m, q, sz_q);
                unlock_ordered(p, q);
            }
            else
            {
                unlock_ordered(p, q);
                // Use find_and_update so that `f` is applied to the existing
                // overflow value rather than blindly overwriting it with `val`.
                const auto r = ov->find_and_update(key,
                    [&f](const Overflow_Map_Val& ov_val) -> Overflow_Map_Val {
                        return {ov_val.cs, ov_val.min_coord, f(ov_val.val)};
                    });
                if(r) return val_t(r->val);
                ov->insert(key, Overflow_Map_Val(c, m, val));
                ov_approx_sz_.fetch_add(1, std::memory_order_relaxed);
                ov_total_inserts_.fetch_add(1, std::memory_order_relaxed);
                tl_ov_happened_ = true;
                return null_val;
            }

            return null_val;
        }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::try_resize()
{
    if constexpr(mt_)   table_lock.lock();

    const uint64_t cur_sz = approx_sz.load(std::memory_order_relaxed);
    const uint64_t cur_ov = ov_approx_sz_.load(std::memory_order_relaxed);
    if(cur_sz >= resize_th || cur_ov >= ov_resize_th_)
    {
        ResizeEvent ev;
        ev.overflow_triggered = (cur_ov >= ov_resize_th_);
        ev.old_cap  = capacity();
        ev.ov_count = cur_ov;
        ev.main_sz  = cur_sz;
        const auto t0 = std::chrono::steady_clock::now();
        resize();
        ev.elapsed_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        resize_log_.push_back(ev);
    }

    if constexpr(mt_)   table_lock.unlock();
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::resize()
{
    const auto old_cap = cap_;
    cap_ *= 2;
    T_new = aligned_allocate<flat_t>(cap_ * B, L1_CACHE_LINE_SIZE);
    M_new = aligned_allocate<Metadata>(cap_, L1_CACHE_LINE_SIZE);
    ov_new = new ov_t(capacity() * of_default);
    std::memset(M_new, 0, cap_ * sizeof(Metadata)); // Checksum `0` denotes absence of a key.

    std::vector<std::thread> W;
    W.reserve(resize_worker_c);
    for(uint64_t i = 0; i < resize_worker_c; ++i)
        W.emplace_back([&](const auto w_id)
        {
            const auto stride = resize_worker_c;
            for(std::size_t b = w_id; b < old_cap; b += stride)
                for(std::size_t slot = 0; slot < B; ++slot)
                    if(M[b].cs[slot])
                        insert_at_resize(T[b * B + slot], M[b].cs[slot], M[b].min_coord[slot]);

            auto it = ov->iterator(resize_worker_c, w_id);
            typedef std::conditional_t<!is_map_, Overflow_Set_Val, Overflow_Map_Val> ov_val_t;
            Kmer<k> key;
            ov_val_t val;
            while(it.next(key, val))
                if constexpr(is_set_)
                    insert_at_resize(key, val.cs, val.min_coord);
                else
                    insert_at_resize(flat_t(key, val.val), val.cs, val.min_coord);
        }, i);

    std::for_each(W.begin(), W.end(), [](auto& w){ w.join(); });


    std::swap(M, M_new), std::swap(T, T_new), std::swap(ov, ov_new);
    deallocate(M_new), deallocate(T_new), delete(ov_new);
    M_new = nullptr, T_new = nullptr, ov_new = nullptr;

    resize_th = capacity() * lf;
    ov_approx_sz_.store(0, std::memory_order_relaxed);
    ov_resize_th_ = static_cast<uint64_t>(capacity() * of_default * 0.5);

}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline auto Streaming_Kmer_Hash_Table<k, mt_, T_, l>::find(const Kmer_Window<k, l>& w) const -> find_ret_t
{
    const auto key = w.v.canonical();
    const auto nt_h = w.minimizer_hash();
    const auto h = XXH3_64bits_withSeed(&nt_h, sizeof(nt_h), kBucketSeed);
    const auto idx_mask = cap_ - 1;

    const auto b = h & idx_mask;
    const auto c = std::max(w.rh.template checksum<8>(), uint64_t(1));    // Avoiding checksum 0 by overloading checksum 1.
    if(const auto r = find(key, c, b))
        return r;

    if(bucket_size(b) < B)
        return null_val;

    // TODO: get the following hashes from the minimizer iterator (rolling hashes), instead of full-blown xxHash computations.
    const auto h_1 = XXH3_64bits_withSeed(&h, sizeof(h), 0);
    const auto h_2 = XXH3_64bits_withSeed(&h, sizeof(h), 1llu << 63);

    const auto p = h_1 & idx_mask;
    const auto q = h_2 & idx_mask;
    if(const auto r_p = find(key, c, p))    return r_p;
    if(const auto r_q = find(key, c, q))    return r_q;

    if(bucket_size(p) == B && bucket_size(q) == B)
        return find_in_overflow(key);

    return null_val;
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
template <typename F>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::for_each(F&& f) const
{
    // Main flat table.
    for(auto it = cbegin(); it != cend(); ++it)
        f(*it);

    // Overflow table.
    typedef std::conditional_t<is_set_, Overflow_Set_Val, Overflow_Map_Val> ov_val_t;
    auto ov_it = ov->iterator();
    Kmer<k> key;
    ov_val_t ov_val;
    while(ov_it.next(key, ov_val))
    {
        if constexpr(is_set_)   f(key);
        else                    { flat_t e(key, ov_val.val); f(e); }
    }
}


template <uint16_t k, bool mt_, typename T_, uint16_t l>
inline void Streaming_Kmer_Hash_Table<k, mt_, T_, l>::report_stats() const
{
    std::cerr   << "Table height: " << cap_ << ".\n";
    std::cerr   << "Table width:  " << B << ".\n";
    std::cerr   << "========================================\n";


    uint64_t sz = size();
    std::cerr   << "Complete table size: " << sz << ".\n";
    std::cerr   << "Overflow table size: " << ov->size() << " ("
                    << (ov->size() * 100.0 / sz) << "%).\n";
    std::cerr   << "========================================\n";


    uint64_t full = 0, empty = 0, blk_occ = 0;
    for(std::size_t b = 0; b < cap_; ++b)
        full += (bucket_size(b) == B), empty += (bucket_size(b) == 0), blk_occ += bucket_size(b);

    std::cerr   << "#Blocks:       " << cap_ << ".\n";
    std::cerr   << "#Full-blocks:  " << full << "(" << (full * 100.0 / cap_) << "%).\n";
    std::cerr   << "#Empty-blocks: " << empty << "(" << (empty * 100.0 / cap_) << "%).\n";
    std::cerr   << "Avg-occupancy: " << blk_occ / cap_ << ".\n";
    std::cerr   << "========================================\n";

}


// Thread-local flag: set by overflow insert path, consumed by public upsert.
template <uint16_t k, bool mt_, typename T_, uint16_t l>
thread_local bool Streaming_Kmer_Hash_Table<k, mt_, T_, l>::tl_ov_happened_ = false;


}



#endif
