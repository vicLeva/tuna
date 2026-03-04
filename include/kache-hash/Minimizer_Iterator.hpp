
#ifndef KACHE_HASH_MINIMIZER_ITERATOR_HPP
#define KACHE_HASH_MINIMIZER_ITERATOR_HPP



#include "DNA.hpp"
#include "DNA_Utility.hpp"
#include "Kmer.hpp"
#include "xxHash/xxhash.h"

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <cassert>


namespace kache_hash
{


// A class to iterate over the `l`-minimizers of the constituent `k`-mers of a
// given sequence, computing the minimizer for each `k`-mer in amortized `O(1)`
// time. Canonical minimizers are computed, i.e. both strand-forms of the
// sequence is considered.
template <uint16_t k, uint16_t l>
class Min_Iterator
{
    static_assert(k <= 64, "k-mer size must not exceed 64 for internal purposes of the hash table.");

private:

    // Returns the minimum of `x` and `y`.
    static constexpr auto min = [](const uint64_t x, const uint64_t y){ return x < y ? x : y; };

    // Returns the 64-bit equivalent of `val`.
    template <typename T_> static constexpr uint64_t as_u64(const T_ val) { return static_cast<uint64_t>(val); }

    static constexpr auto cmov = [](const bool cond, const int16_t x, const int16_t y) { return cond ? x : y; };

    static constexpr std::size_t fwd = 0;   // Marker for data in the forward-strand.
    static constexpr std::size_t rev = 1;   // Marker for data in the reverse-strand.

    static constexpr uint16_t w = k - l + 1;    // `l`-mer window width in a `k`-mer, i.e. number of `l`-mers in a `k`-mer.

    uint64_t last_lmer[2];  // The last l-mer processed from the sequence, in the forward- and in the reverse-strand.
    static constexpr uint64_t clear_MSN_mask = ~(as_u64(0b11) << (2 * (l - 1)));    // Bitmask to clear the most-significant nucleotide bits of l-mers.

    uint64_t H[2][w + 1];   // The hashes of the l-mers in a k-mer window. The l-mers sit at `H[0..k - l]`, and is conceptually broken into two sub-windows:
                            // a prefix and a suffix one.

    uint64_t M_pre[2][w + 1];   // Suffix-mins of the prefix window.
    uint64_t M_suf[2];  // Min of the suffix window.

    uint16_t P_pre[2][w + 1];   // Positions (absolute) of the suffix-mins at the prefix window.
    uint16_t P_suf[2];  // Position (absolute) of the min at the suffix window.

    std::size_t pivot;  // Pivot-index, inclusive to the prefix window, breaking the window into the two sub-windows.

    // Moves the suffix window to the prefix window, and empties it—resetting
    // them to their default condition. Also recomputes the aggregate minimum
    // results of the windows.
    void reset_windows();


public:

    // Constructs a minimizer iterator to iterate over `l`-minimizers of
    // the `k`-mers of a given sequence in a streaming manner.
    Min_Iterator() {}

    Min_Iterator(const Min_Iterator&) = delete;
    Min_Iterator(Min_Iterator&&) = delete;
    Min_Iterator& operator=(const Min_Iterator&) = delete;
    Min_Iterator& operator=(Min_Iterator&&) = delete;

    // Returns the mask to extract the strand-orientation of a k-mer from
    // its minimizer-coordinate, denoting where the minimizer comes from.
    static constexpr uint8_t min_orientation_mask() { return 0b0100'0000; }

    // Returns the mask to extract the index of the minimizer from the
    // appropriate strand of the k-mer.
    static constexpr uint8_t min_index_mask() { return 0b0011'1111; }

    // Returns the 64-bit hash value of the `l`-mer `lmer`.
    static uint64_t lmer_hash(uint64_t lmer);

    // Resets the iterator to the sequence `seq`. The iterator sits at the
    // first `k`-mer after the construction.
    void reset(const char* seq);

    // Resets the iterator to the k-mer `kmer`.
    void reset(const Kmer<k>& kmer);

    // Returns the minimizer's 64-bit hash value.
    uint64_t hash() const;

    // Returns the minimizer's 64-bit hash value. Sets the second-highest bit of
    // `min_coord` to `1` iff the minimizer comes from the forward strand. Also
    // places the position of the minimizer `l`-mer in its corresponding strand
    // at `min_coord`'s lowest 6 bits.
    uint64_t hash(uint8_t& min_coord) const;

    // Moves the iterator to a next k-mer extending it with the nucleobase `b`.
    void advance(DNA::Base b);

    // Moves the iterator to a next k-mer extending it with the character `ch`.
    void advance(char ch);

    // Returns the minimizer's position in the `k`-mer in the forward-strand.
    // If the minimizer actually comes from the reverse-strand, then the l-mer
    // at this index must be taken in the RC-orientation to get the minimizer.
    // If there are multiple candidates, then returns some value `< 0`.
    // TODO: rename.
    int16_t hash_pos() const;

    // Returns the minimizer-hash of the k-mer `kmer`.
    static uint64_t hash(Kmer<k> kmer);
};


template <uint16_t k, uint16_t l>
inline uint64_t Min_Iterator<k, l>::lmer_hash(const uint64_t lmer)
{
    const auto h = [&](const uint64_t lmer)
        {
            // TODO: does making the size of the key unknown at compile-time non-trivially affect runtime?
            // Can use template-metaprogramming to solve this.
            // return XXH3_64bits_withSeed(&lmer, lmer_byte_c, min_seed);
            constexpr uint64_t min_seed = 0;
            return XXH3_64bits_withSeed(&lmer, sizeof(lmer), min_seed); // TODO: replace with rolling hash: `xxHash` is too costly in this context.
        };

    return h(lmer);
}


template <uint16_t k, uint16_t l>
inline void Min_Iterator<k, l>::reset(const char* const seq)
{
    last_lmer[fwd] = last_lmer[rev] = 0;

    for(std::size_t idx = 0; idx < l; ++idx)
    {
        assert(DNA_Utility::is_DNA_base(seq[idx]));

        const auto b = DNA_Utility::map_base(seq[idx]);
        assert(b != DNA::N);
        last_lmer[fwd] |= (as_u64(b) << (2 * (l - 1 - idx)));
        last_lmer[rev] |= (as_u64(DNA_Utility::complement(b)) << (2 * idx));
    }

    pivot = 0;
    H[fwd][pivot] = lmer_hash(last_lmer[fwd]);
    H[rev][pivot] = lmer_hash(last_lmer[rev]);

    for(std::size_t idx = l; idx < k; ++idx)
    {
        assert(DNA_Utility::is_DNA_base(seq[idx]));

        const DNA::Base b = DNA_Utility::map_base(seq[idx]);
        assert(b != DNA::N);
        last_lmer[fwd] = ((last_lmer[fwd] & clear_MSN_mask) << 2) | b;
        last_lmer[rev] = (last_lmer[rev] >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

        pivot++;
        H[fwd][pivot] = lmer_hash(last_lmer[fwd]);
        H[rev][pivot] = lmer_hash(last_lmer[rev]);
    }

    reset_windows();
}


template <uint16_t k, uint16_t l>
inline void Min_Iterator<k, l>::reset(const Kmer<k>& kmer)
{
    last_lmer[fwd] = last_lmer[rev] = 0;

    for(std::size_t idx = 0; idx < l; ++idx)
    {
        const auto b = DNA_Utility::map_base(kmer.base_at(k - 1 - idx));
        last_lmer[fwd] |= (as_u64(b) << (2 * (l - 1 - idx)));
        last_lmer[rev] |= (as_u64(DNA_Utility::complement(b)) << (2 * idx));
    }

    pivot = 0;
    H[fwd][pivot] = lmer_hash(last_lmer[fwd]);
    H[rev][pivot] = lmer_hash(last_lmer[rev]);

    for(std::size_t idx = l; idx < k; ++idx)
    {
        const DNA::Base b = DNA_Utility::map_base(kmer.base_at(k - 1 - idx));
        last_lmer[fwd] = ((last_lmer[fwd] & clear_MSN_mask) << 2) | b;
        last_lmer[rev] = (last_lmer[rev] >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

        pivot++;
        H[fwd][pivot] = lmer_hash(last_lmer[fwd]);
        H[rev][pivot] = lmer_hash(last_lmer[rev]);
    }

    reset_windows();
}


template <uint16_t k, uint16_t l>
inline void Min_Iterator<k, l>::reset_windows()
{
    assert(pivot == w - 1);

    // TODO: try vectorization, if not already done by the compiler.

    M_pre[fwd][w] = M_pre[rev][w] = std::numeric_limits<uint64_t>::max();
    P_pre[fwd][w] = P_pre[rev][w] = w;
    for(int16_t i = w - 1; i >= 0; --i)
    {
        M_pre[fwd][i] = min(M_pre[fwd][i + 1], H[fwd][i]),
        P_pre[fwd][i] = cmov(M_pre[fwd][i + 1] <= H[fwd][i], P_pre[fwd][i + 1], i);

        M_pre[rev][i] = min(M_pre[rev][i + 1], H[rev][i]),
        P_pre[rev][i] = cmov(M_pre[rev][i + 1] <= H[rev][i], P_pre[rev][i + 1], i);
    }

    M_suf[fwd] = M_suf[rev] = std::numeric_limits<uint64_t>::max();
    P_suf[fwd] = P_suf[rev] = w;
    pivot = 0;
}


template <uint16_t k, uint16_t l>
inline uint64_t Min_Iterator<k, l>::hash() const
{
    return min(min(M_pre[fwd][pivot], M_suf[fwd]), min(M_pre[rev][pivot], M_suf[rev]));
}


template <uint16_t k, uint16_t l>
inline uint64_t Min_Iterator<k, l>::hash(uint8_t& min_coord) const
{
    const auto h_fwd = min(M_pre[fwd][pivot], M_suf[fwd]);
    const auto h_rev = min(M_pre[rev][pivot], M_suf[rev]);

    if(h_fwd <= h_rev)
    {
        min_coord = cmov(M_pre[fwd][pivot] < M_suf[fwd], P_pre[fwd][pivot] - pivot, w - pivot + P_suf[fwd]) | min_orientation_mask();
        return h_fwd;
    }
    else
    {
        min_coord = k - 1 - (cmov(M_pre[rev][pivot] < M_suf[rev], P_pre[rev][pivot] - pivot, w - pivot + P_suf[rev]) + l - 1);
        return h_rev;
    }
}


template <uint16_t k, uint16_t l>
inline int16_t Min_Iterator<k, l>::hash_pos() const
{
    const auto h_fwd = min(M_pre[fwd][pivot], M_suf[fwd]);
    const auto h_rev = min(M_pre[rev][pivot], M_suf[rev]);
    return  h_fwd <= h_rev ?
                cmov(M_pre[fwd][pivot] < M_suf[fwd],
                    P_pre[fwd][pivot] - pivot,
                    w - pivot + P_suf[fwd]) :
                cmov(M_pre[rev][pivot] < M_suf[rev],
                    P_pre[rev][pivot] - pivot,
                    w - pivot + P_suf[rev]) + l - 1;
}


template <uint16_t k, uint16_t l>
inline void Min_Iterator<k, l>::advance(const DNA::Base b)
{
    assert(b != DNA::N);
    last_lmer[fwd] = ((last_lmer[fwd] & clear_MSN_mask) << 2) | b;
    last_lmer[rev] = (last_lmer[rev] >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

    H[fwd][pivot] = lmer_hash(last_lmer[fwd]);
    H[rev][pivot] = lmer_hash(last_lmer[rev]);

    P_suf[fwd] = cmov(M_suf[fwd] < H[fwd][pivot], P_suf[fwd], pivot);
    M_suf[fwd] = min(M_suf[fwd], H[fwd][pivot]);

    P_suf[rev] = cmov(M_suf[rev] < H[rev][pivot], P_suf[rev], pivot);
    M_suf[rev] = min(M_suf[rev], H[rev][pivot]);

    if(pivot < w - 1)
        pivot++;
    else
        reset_windows();
}


template <uint16_t k, uint16_t l>
inline void Min_Iterator<k, l>::advance(const char ch)
{
    assert(DNA_Utility::is_DNA_base(ch));
    advance(DNA_Utility::map_base(ch));
}


template <uint16_t k, uint16_t l>
inline uint64_t Min_Iterator<k, l>::hash(const Kmer<k> kmer)
{
    uint64_t last_lmer[2] = {0, 0};

    for(std::size_t idx = 0; idx < l; ++idx)
    {
        const auto b = kmer.base_at(k - 1 - idx);
        last_lmer[fwd] |= (as_u64(b) << (2 * (l - 1 - idx)));
        last_lmer[rev] |= (as_u64(DNA_Utility::complement(b)) << (2 * idx));
    }

    auto h_fwd = lmer_hash(last_lmer[fwd]);
    auto h_rev = lmer_hash(last_lmer[rev]);
    for(std::size_t idx = l; idx < k; ++idx)
    {
        const DNA::Base b = kmer.base_at(k - 1 - idx);
        last_lmer[fwd] = ((last_lmer[fwd] & clear_MSN_mask) << 2) | b;
        last_lmer[rev] = (last_lmer[rev] >> 2) | (as_u64(DNA_Utility::complement(b)) << (2 * (l - 1)));

        h_fwd = min(h_fwd, lmer_hash(last_lmer[fwd]));
        h_rev = min(h_rev, lmer_hash(last_lmer[rev]));
    }

    return min(h_fwd, h_rev);
}


}



#endif
