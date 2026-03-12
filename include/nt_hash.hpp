#pragma once

// Canonical ntHash rolling hash for DNA l-mers.
//
// Reference: Mohamadi et al., "ntHash: recursive nucleotide hashing",
//            Bioinformatics 2016. Classic ntHash (not ntHash2).
//
// Canonical hash = fwd_hash XOR rev_hash, matching the convention used by
// simd-minimizers (Rust crate seq-hash, NtHasher<true>).
//
// Usage:
//   nt_hash::Roller h(l);
//   h.init(seq);                       // initialise on seq[0..l-1]
//   uint64_t hc = h.canonical();
//   h.roll(out_2bit, in_2bit);         // slide window by 1
//   uint64_t hc2 = h.canonical();

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

inline uint64_t rol64(uint64_t x, int s) noexcept { return (x << s) | (x >> (64 - s)); }
inline uint64_t ror64(uint64_t x, int s) noexcept { return (x >> s) | (x << (64 - s)); }


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
// Maintains the ntHash of a sliding l-mer window.
//
// Forward hash formula:
//   H_fwd(s[0..l-1]) = XOR_{i=0}^{l-1} rol(FWD[s[i]], l-1-i)
// Rolling update (drop s_out from left, add s_in on right):
//   H_fwd_new = rol(H_fwd, 1) ^ rol(FWD[s_out], l) ^ FWD[s_in]
//
// Reverse-complement hash formula:
//   H_rev(s[0..l-1]) = XOR_{i=0}^{l-1} rol(REV[s[i]], i)
//   (this equals the forward ntHash of the RC l-mer comp(s[l-1])…comp(s[0]))
// Rolling update:
//   H_rev_new = ror(H_rev, 1) ^ ror(REV[s_out], 1) ^ rol(REV[s_in], l-1)

class Roller {
    uint64_t fwd_ = 0, rev_ = 0;
    uint16_t l_;

public:
    explicit Roller(uint16_t l) noexcept : l_(l) {}

    // Initialise from seq[0..l-1].
    void init(const char* seq) noexcept {
        fwd_ = 0; rev_ = 0;
        for (uint16_t i = 0; i < l_; ++i) {
            const uint8_t b = to_2bit(seq[i]);
            fwd_ ^= rol64(FWD[b], l_ - 1 - i);
            rev_ ^= rol64(REV[b], i);
        }
    }

    // Slide the window: `out_2bit` leaves from the left, `in_2bit` enters on the right.
    void roll(uint8_t out_2bit, uint8_t in_2bit) noexcept {
        fwd_ = rol64(fwd_, 1) ^ rol64(FWD[out_2bit], l_) ^ FWD[in_2bit];
        rev_ = ror64(rev_, 1) ^ ror64(REV[out_2bit], 1) ^ rol64(REV[in_2bit], l_ - 1);
    }

    uint64_t fwd()       const noexcept { return fwd_; }
    uint64_t rev()       const noexcept { return rev_; }
    // Canonical hash = fwd XOR rev (simd-minimizers convention).
    uint64_t canonical() const noexcept { return fwd_ ^ rev_; }
};

} // namespace nt_hash
