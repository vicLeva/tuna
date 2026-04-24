#pragma once

// Forward-only ntHash-style rolling hash for protein m-mers.
// No reverse complement (proteins are directional).
// 20 per-amino-acid random seed constants; O(1) amortised per advance.

#include "aa_encode.hpp"
#include <array>
#include <cstdint>

namespace prot {

// 20 pseudo-random 64-bit seeds, one per amino acid (A…Y, alphabetical order).
// Generated from known good bit-mixing constants.
static constexpr uint64_t AA_SEEDS[20] = {
    0x3c8bfbb395c60474ULL, 0x3193c18562a02b4cULL, 0x295549f54be24456ULL,
    0x20323ed082572324ULL, 0x6c62272e07bb0142ULL, 0x62b821756295c58dULL,
    0x516e3d1c1c9d0e89ULL, 0x4d4c6e7ab3fa18e2ULL, 0x38b8b0b2f6c92c5dULL,
    0xa09e667f3bcc908bULL, 0xb452821e638d0135ULL, 0xbe67572b08a8e94bULL,
    0x985f24a9638c2d5eULL, 0xc89ca6d6584a0fceULL, 0xd5a60faa61b2e34eULL,
    0xe14e5b5c46de7f04ULL, 0xf2c56c3eb7e08b29ULL, 0x17a58b0b84fd11d8ULL,
    0x279dba57b6d51a39ULL, 0x2a9bc8b8b5b3a7c6ULL,
};

inline constexpr uint64_t aa_rol64(uint64_t x, int s) noexcept {
    return (x << s) | (x >> (64 - s));
}

// Rolling hash for a window of m amino acids.
// H = XOR_i rol(SEED[enc[i]], m-1-i)
// roll: H = rol(H, 1) ^ rol(SEED[out], m) ^ SEED[in]
template <uint16_t m>
class ProtRoller {
    static_assert(m >= 2 && m <= 15, "protein minimizer m must be in [2, 15]");

    // CTAB[out * 20 + in] = rol(SEED[out], m) ^ SEED[in]
    // 400 entries × 8 B = 3.2 KB — stays in L1 cache.
    static constexpr std::array<uint64_t, 400> make_ctab() noexcept {
        std::array<uint64_t, 400> t{};
        for (int out = 0; out < 20; ++out)
            for (int in = 0; in < 20; ++in)
                t[out * 20 + in] = aa_rol64(AA_SEEDS[out], m) ^ AA_SEEDS[in];
        return t;
    }
    static constexpr auto CTAB = make_ctab();

    uint64_t h_ = 0;

public:
    ProtRoller() noexcept = default;

    // Initialise from m encoded amino acids (5-bit values 0-19).
    void init(const uint8_t* enc) noexcept {
        h_ = 0;
        for (uint16_t i = 0; i < m; ++i)
            h_ ^= aa_rol64(AA_SEEDS[enc[i]], m - 1 - i);
    }

    // Slide window: out5 leaves on the left, in5 enters on the right.
    void roll(uint8_t out5, uint8_t in5) noexcept {
        h_ = aa_rol64(h_, 1) ^ CTAB[out5 * 20 + in5];
    }

    uint64_t hash() const noexcept { return h_; }
};

} // namespace prot
