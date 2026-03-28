#pragma once

// Canonical ntHash rolling hash for DNA m-mers.
//
// Reference: Mohamadi et al., "ntHash: recursive nucleotide hashing",
//            Bioinformatics 2016. Classic ntHash (not ntHash2).
//
// Canonical hash = fwd_hash XOR rev_hash, matching the convention used by
// simd-minimizers (Rust crate seq-hash, NtHasher<true>).
//
// Usage:
//   nt_hash::Roller h(m);
//   h.init(seq);                       // initialise on seq[0..m-1]
//   uint64_t hc = h.canonical();
//   h.roll(out_2bit, in_2bit);         // slide window by 1
//   uint64_t hc2 = h.canonical();

#include <array>
#include <cstdint>

namespace nt_hash {

// ── Per-nucleotide seeds (A=0, C=1, T=2, G=3) ────────────────────────────────
// Taken from the ntHash reference implementation.
// Encoding matches the raw bit trick (c >> 1) & 3 on ASCII ACGT.
static constexpr uint64_t FWD[4] = {
    0x3c8bfbb395c60474ULL, // A
    0x3193c18562a02b4cULL, // C
    0x295549f54be24456ULL, // T
    0x20323ed082572324ULL, // G
};

// Reverse-complement seeds: REV[b] = FWD[complement(b)], where A(0)↔T(2), C(1)↔G(3).
// complement = b ^ 2.
static constexpr uint64_t REV[4] = {
    0x295549f54be24456ULL, // complement(A) = T
    0x20323ed082572324ULL, // complement(C) = G
    0x3c8bfbb395c60474ULL, // complement(T) = A
    0x3193c18562a02b4cULL, // complement(G) = C
};


// ── Bit helpers ───────────────────────────────────────────────────────────────

constexpr uint64_t rol64(uint64_t x, int s) noexcept { return (x << s) | (x >> (64 - s)); }
constexpr uint64_t ror64(uint64_t x, int s) noexcept { return (x >> s) | (x << (64 - s)); }


// ── DNA helpers ───────────────────────────────────────────────────────────────

// Map ASCII DNA (ACGT, upper or lower case) to 2-bit (A=0, C=1, T=2, G=3).
// Behaviour is undefined for non-ACGT input.
inline uint8_t to_2bit(char c) noexcept {
    return (static_cast<uint8_t>(c) >> 1) & 3u;
}

// 2-bit complement: A(0)↔T(2), C(1)↔G(3).
inline uint8_t complement_2bit(uint8_t b) noexcept { return b ^ 2u; }

// Returns true iff c is one of ACGT (upper or lower case).
// Used to detect placeholder/ambiguous characters (N, gaps, …).
inline bool is_dna(char c) noexcept {
    const uint8_t u = static_cast<uint8_t>(c) | 0x20u; // fold to lowercase
    return u == 'a' || u == 'c' || u == 'g' || u == 't';
}


// ── Rolling hash ─────────────────────────────────────────────────────────────
//
// Maintains forward and reverse-complement ntHash of a sliding m-mer window.
//   H_fwd = XOR_{i} rol(FWD[s[i]], m-1-i); roll: rol(H_fwd,1) ^ rol(FWD[out],m) ^ FWD[in]
//   H_rev = XOR_{i} rol(REV[s[i]], i);      roll: ror(H_rev,1) ^ ror(REV[out],1) ^ rol(REV[in],m-1)

template <uint16_t m>
class Roller {
    // Precomputed roll contributions indexed by (out_2bit << 2) | in_2bit.
    // CTAB_FWD[idx] = rol(FWD[out], m) ^ FWD[in]
    // CTAB_REV[idx] = ror(REV[out], 1) ^ rol(REV[in], m-1)
    // Halves the number of table lookups/rotations in roll(), keeping fwd_/rev_ in GPRs.
    static constexpr std::array<uint64_t, 16> make_ctab_fwd() noexcept {
        std::array<uint64_t, 16> t{};
        for (int out = 0; out < 4; ++out)
            for (int in  = 0; in  < 4; ++in )
                t[out * 4 + in] = rol64(FWD[out], m) ^ FWD[in];
        return t;
    }
    static constexpr std::array<uint64_t, 16> make_ctab_rev() noexcept {
        std::array<uint64_t, 16> t{};
        for (int out = 0; out < 4; ++out)
            for (int in  = 0; in  < 4; ++in )
                t[out * 4 + in] = ror64(REV[out], 1) ^ rol64(REV[in], m - 1);
        return t;
    }

    static constexpr auto CTAB_FWD = make_ctab_fwd();
    static constexpr auto CTAB_REV = make_ctab_rev();

    uint64_t fwd_ = 0, rev_ = 0;

public:
    Roller() noexcept = default;

    // Initialise from seq[0..m-1].
    void init(const char* seq) noexcept {
        fwd_ = 0; rev_ = 0;
        for (uint16_t i = 0; i < m; ++i) {
            const uint8_t b = to_2bit(seq[i]);
            fwd_ ^= rol64(FWD[b], m - 1 - i);
            rev_ ^= rol64(REV[b], i);
        }
    }

    // Slide the window: `out_2bit` leaves from the left, `in_2bit` enters on the right.
    void roll(uint8_t out_2bit, uint8_t in_2bit) noexcept {
        const uint8_t idx = static_cast<uint8_t>((out_2bit << 2) | in_2bit);
        fwd_ = rol64(fwd_, 1) ^ CTAB_FWD[idx];
        rev_ = ror64(rev_, 1) ^ CTAB_REV[idx];
    }

    uint64_t fwd()       const noexcept { return fwd_; }
    uint64_t rev()       const noexcept { return rev_; }
    // Canonical hash = fwd XOR rev (simd-minimizers convention).
    uint64_t canonical() const noexcept { return fwd_ ^ rev_; }
};

} // namespace nt_hash
