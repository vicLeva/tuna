
#ifndef KACHE_HASH_ROLLING_HASH_HPP
#define KACHE_HASH_ROLLING_HASH_HPP



#include "DNA.hpp"
#include "DNA_Utility.hpp"
#include "Kmer.hpp"

#include <cstdint>
#include <cstddef>
#include <cassert>


namespace kache_hash
{


// Rolling tabulation hasher for `k`-mers, using cyclic polynomials.
template <uint16_t k, bool canonical = false>
class Rolling_Hash
{
private:

    uint8_t base[k];    // Circular queue of the nucleobases in the current k-mer.
    std::size_t off;    // Index of the front of the queue.

#ifdef DEVELOP
    static constexpr auto t = "12:34:56";
#else
    static constexpr auto t = __TIME__;
#endif
    static constexpr auto atoi = [](const char ch){ return static_cast<uint64_t>(ch - '0'); };

    static constexpr uint64_t a = 6364136223846793005llu;   // Multiplier in LCG: from Knuth's MMIX.
    static constexpr uint64_t c = 1442695040888963407llu;   // Increment in LCG: from Knuth's MMIX.
    static constexpr uint64_t s_a = (atoi(t[0]) * 10 + atoi(t[1])) * 3600 + (atoi(t[3]) * 10 + atoi(t[4])) * 60 + (atoi(t[6]) * 10 + atoi(t[7]));
    static constexpr uint64_t s_c = a * s_a + c;
    static constexpr uint64_t s_g = a * s_c + c;
    static constexpr uint64_t s_t = a * s_g + c;
    static constexpr uint64_t s[] = {s_a, s_c, s_g, s_t};   // Substitution function for nucleobases.

    // Returns the rotation of `x` to the left by `r` bits.
    static constexpr uint64_t rotl(const uint64_t x, const uint32_t r) { return (r == 0) ? x : (x << r) | (x >> (64 - r)); }

    // Returns the rotation of `x` to the right by `r` bits.
    static constexpr uint64_t rotr(const uint64_t x, const uint32_t r) { return (x >> r) | (x << (64 - r)); }

    static constexpr uint64_t rotl_k[]   = {rotl(s_a, k), rotl(s_c, k), rotl(s_g, k), rotl(s_t, k)};    // Left-rotation result by `k`.
    static constexpr uint64_t rotl_km1[] = {rotl(s_a, k - 1), rotl(s_c, k - 1), rotl(s_g, k - 1), rotl(s_t, k - 1)};    // Left-rotation result by `k - 1`.
    static constexpr uint64_t rotr_1[]   = {rotr(s_a, 1), rotr(s_c, 1), rotr(s_g, 1), rotr(s_t, 1)};    // Right-rotation result.

    uint64_t h_f;   // Hash of the current k-mer in the forward-strand.
    uint64_t h_r;   // Hash of the current k-mer in the reverse-strand.

public:

    // Constructs an empty rolling tabulation hasher. This must be initialized
    // before use.
    Rolling_Hash() {}

    // Constructs a rolling tabulation hasher for `k`-mers initialized at the
    // beginning of the sequence `seq`.
    Rolling_Hash(const char* seq) { init(seq); }

    // Initializes the hasher for `k`-mers initialized at the beginning of the
    // sequence `seq`.
    void init(const char* seq);

    // Initializes the hasher for the `k`-mer `kmer`. The minimizer-hash is not
    // initialized and is assumed to be computed externally.
    void init(const Kmer<k>& kmer);

    // Advances the hasher in its underlying sequence by the nucleobase `b` to
    // the right.
    void advance(DNA::Base b);

    // Advances the hasher in its underlying sequence by the character `ch` to
    // the right.
    void advance(char ch);

    // Returns the hash of the current k-mer in the forward-strand.
    auto hash_fwd() const { return h_f; }

    // Returns the hash of the current k-mer in the reverse-strand.
    auto hash_rev() const { return h_r; }

    // Returns the current hash.
    uint64_t hash() const;

    // Returns an `N`-bit checksum of the hash.
    template <uint8_t N = 8> uint64_t checksum() const { return hash() >> (64 - N); }
};


template <uint16_t k, bool canonical>
inline void Rolling_Hash<k, canonical>::init(const char* const seq)
{
    h_f = h_r = 0;
    off = 0;
    for(std::size_t i = 0; i < k; ++i)
    {
        base[i] = DNA_Utility::map_base(seq[i]);
        h_f ^= rotl(s[base[i]], k - 1 - i);
        if constexpr(canonical)
            h_r ^= rotl(s[DNA_Utility::complement(DNA::Base(base[i]))], i);
    }
}


template <uint16_t k, bool canonical>
inline void Rolling_Hash<k, canonical>::init(const Kmer<k>& kmer)
{
    h_f = h_r = 0;
    off = 0;
    for(std::size_t i = 0; i < k; ++i)
    {
        base[i] = kmer.base_at(k - 1 - i);
        h_f ^= rotl(s[base[i]], k - 1 - i);
        if constexpr(canonical)
            h_r ^= rotl(s[DNA_Utility::complement(DNA::Base(base[i]))], i);
    }
}


template <uint16_t k, bool canonical>
inline void Rolling_Hash<k, canonical>::advance(const DNA::Base b)
{
    const auto out = base[off];
    base[off] = b;
    const auto in = base[off];
    off = (off + 1 >= k ? 0 : off + 1);

    h_f = rotl(h_f ^ rotl_km1[out], 1) ^ s[in];
    if constexpr(canonical)
    {
        const auto o = DNA_Utility::complement(DNA::Base(out));
        const auto i = DNA_Utility::complement(DNA::Base(in));
        h_r = rotr(h_r ^ s[o], 1) ^ rotl_km1[i];
    }
}


template <uint16_t k, bool canonical>
inline void Rolling_Hash<k, canonical>::advance(const char ch)
{
    advance(DNA_Utility::map_base(ch));
}


template <uint16_t k, bool canonical>
inline uint64_t Rolling_Hash<k, canonical>::hash() const
{
    if constexpr(canonical)
        return h_f + h_r;

    return h_f;
}


}



#endif
