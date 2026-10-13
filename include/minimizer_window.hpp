#pragma once

// MinimizerWindow<k, m> — streaming m-minimizer iterator, sshash scheme.
//
// The order over m-mers is the one sshash uses (include/minimizer_iterator.hpp,
// util::canonical_mmer_at): the canonical m-mer is the smaller of the m-mer and
// its reverse complement as packed 2-bit values, and that value is hashed. This
// is a random minimizer scheme over the strand-symmetric classes, so the usual
// density results apply to it.
//
// The previous order hashed each strand with ntHash and XORed the two. That is
// also strand-symmetric and was cheaper, but it is not "hash of the canonical
// m-mer", so the published density of the random minimizer scheme did not
// describe it.
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over those values.
// Window size w = k - m; arrays sized exactly w+1. O(1) amortised per advance().
//
// Usage:
//   MinimizerWindow<k, m> win;
//   win.reset(seq);          // init from first k ASCII characters
//   win.advance(ch);         // slide by one character
//   win.hash();              // hash of the current window's canonical minimizer

#include "nt_hash.hpp"

#include <cstdint>
#include <limits>

// Partition routing needs a well-mixed value, and the minimizer hash is not
// one. sshash's mixer is a single multiply, so its low bits depend only on the
// low bits of the m-mer, while the high bits of a window *minimum* are pushed
// toward zero by the minimum itself. Measured on one E. coli genome into 4096
// partitions, routing straight off the minimizer gives chi2/df = 114 on the low
// bits and 7685 on the high bits (2092 partitions empty); through this
// finalizer it is 1.18. It runs once per superkmer, not once per base.
inline uint64_t minimizer_mix(uint64_t x) noexcept {
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27; x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    return x;
}

template <uint16_t k, uint16_t m>
class MinimizerWindow {
    static_assert(m < k,  "minimizer length m must be strictly less than k-mer length k");
    // lmer_ stores the 2-bit packed m-mer in a uint64_t — limits m to 32.
    static_assert(m <= 32, "minimizer length m must be <= 32 (lmer_ is uint64_t; 2*m bits required)");

    // Window size: number of m-mers in a k-mer window.
    static constexpr uint16_t w = k - m;

    // Mask to clear the top 2 bits of the 2*m-bit m-mer (used in advance_impl).
    static constexpr uint64_t clear_msn_ = ~(uint64_t(0b11) << (2 * (m - 1)));

    // 2-bit packed m-mer in nt encoding (A=0,C=1,T=2,G=3), forward and reverse
    // complement, kept in step. In this encoding complementing a base is b^2
    // (A=0<->T=2, C=1<->G=3), so the reverse complement rolls in two operations:
    // shift right by one base and bring the complement of the incoming base in
    // at the top.
    uint64_t lmer_ = 0;
    uint64_t rc_lmer_ = 0;

    // sshash's mixer_64 (include/hash_util.hpp): one multiply, one xor. The
    // magic is a seed and is not optional here: without it a poly-A m-mer, which
    // packs to 0, would hash to 0 and so win every window it appears in.
    static constexpr uint64_t MIX_MUL = 0x517cc1b727220a95ull;
    static constexpr uint64_t MIX_MAGIC = 0x9e3779b97f4a7c15ull;

    // Hash of the canonical m-mer currently in lmer_/rc_lmer_.
    uint64_t mmer_hash() const noexcept {
        const uint64_t canon = lmer_ < rc_lmer_ ? lmer_ : rc_lmer_;
        return (canon * MIX_MUL) ^ MIX_MAGIC;
    }

    // Roll both strands by one nt-encoded base.
    void roll_mmer(uint8_t in_2bit) noexcept {
        lmer_    = ((lmer_ & clear_msn_) << 2) | in_2bit;
        rc_lmer_ = (rc_lmer_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (m - 1)));
    }


    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w]: circular buffer of w+1 hash values. pivot = next write index
    //   (= leftmost m-mer, evicted on the next advance).
    // M_pre[i] = min(H[0..i])  — prefix minimums from the rightmost end.
    // M_suf    = running min of hashes written since the last reset_windows().
    // hash() = min(M_pre[pivot], M_suf).  pivot reaches 0 → reset_windows().
    //
    // Minimizer positions are deliberately not tracked: superkmers are cut on
    // hash inequality alone, so the index bookkeeping they need was never read.

    uint64_t    H[w + 1];
    uint64_t    M_pre[w + 1];
    uint64_t    M_suf = 0;
    std::size_t pivot = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    void reset_windows() noexcept {
        M_pre[0] = H[0];
        for (std::size_t i = 1; i <= w; ++i)
            // Branchless CMOV: 50%-taken comparisons cause mispredictions otherwise.
            M_pre[i] = umin(H[i], M_pre[i - 1]);
        M_suf = U64_MAX;
        pivot = w;
    }

public:

    MinimizerWindow() noexcept = default;

    // Initialise to the first k characters of `seq` (ASCII ACGT, any case).
    // `seq` must contain at least k valid DNA characters.
    // After reset(), positions are 0-indexed from seq[0].
    void reset(const char* seq) noexcept {
        lmer_ = 0; rc_lmer_ = 0;
        for (uint16_t i = 0; i < m; ++i) {
            const uint8_t b = nt_hash::to_2bit(seq[i]);
            lmer_    |= uint64_t(b)        << (2 * (m - 1 - i));
            rc_lmer_ |= uint64_t(b ^ 2u)   << (2 * i);
        }

        pivot    = w;
        H[pivot] = mmer_hash();         // m-mer seq[0..m-1] at position 0

        for (uint16_t i = m; i < k; ++i) {
            roll_mmer(nt_hash::to_2bit(seq[i]));
            --pivot;
            H[pivot] = mmer_hash();     // m-mer seq[i-m+1..i] at position i-m+1
        }
        // After loop: pivot=0, H[0] = hash of the last m-mer in the initial k-mer window.
        reset_windows();
    }

    // Slide by one ASCII character.
    void advance(char ch) noexcept {
        roll_mmer(nt_hash::to_2bit(ch));

        const uint64_t h = mmer_hash();
        H[pivot] = h;
        M_suf    = umin(h, M_suf);

        if (pivot > 0)
            --pivot;
        else
            reset_windows();
    }

    // Hash of the canonical m-minimizer of the current k-mer window.
    uint64_t hash() const noexcept {
        return umin(M_pre[pivot], M_suf);
    }

    // Same as hash(), but also sets min_coord = 0 (position/orientation not
    // tracked; provided for compatibility with the kache-hash Kmer_Window API).
    uint64_t hash(uint8_t& min_coord) const noexcept {
        min_coord = 0;
        return hash();
    }
};
