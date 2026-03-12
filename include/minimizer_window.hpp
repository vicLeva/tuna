#pragma once

// MinimizerWindow<k> — streaming l-minimizer iterator using canonical ntHash.
//
// Used in:
//   Phase 0+1 (partition_hash.hpp, partition_kmtricks.hpp) — ASCII input path.
//   Phase 2   (kache-hash Kmer_Window) — field replaces Min_Iterator<k,l> (XXH3),
//             but bypassed for hash/kmtricks superkmers (min_pos byte used instead);
//             only active for KMC-mode superkmers (min_pos == 0xFF fallback).
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ntHash values.
// O(k-l) active entries (arrays sized k for compile-time template sizing),
// O(1) amortised per advance() call.
//
// Two input encodings are supported:
//   ASCII  — advance(char ch)         uses nt_hash::to_2bit() (A=0,C=1,T=2,G=3)
//   kache  — advance_kache(uint8_t b) accepts kache encoding (A=0,C=1,G=2,T=3),
//            converts to nt via b ^ (b >> 1), then calls advance_impl().
//
// Public interface:
//   MinimizerWindow<k> win(l);
//   win.reset(seq);              // initialise to first k ASCII characters
//   uint64_t h = win.hash();     // canonical ntHash of the l-minimizer
//   win.advance(ch);             // slide by one ASCII character
//   win.advance_kache(b);        // slide by one kache-encoded base
//   win.hash(min_coord);         // hash + dummy min_coord (for kache-hash API)

#include "nt_hash.hpp"

#include <cstdint>
#include <cassert>
#include <limits>

template <uint16_t k>
class MinimizerWindow {

    uint16_t l_;
    const uint64_t clear_msn_;  // mask to clear the top 2 bits of the 2*l-bit l-mer

    // 2-bit packed l-mer (forward and RC), tracked to extract the outgoing
    // base cheaply when rolling the ntHash.  nt encoding (A=0,C=1,T=2,G=3).
    uint64_t lmer_     = 0;
    uint64_t lmer_bar_ = 0;

    nt_hash::Roller hasher_;

    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w] is a circular overwrite buffer, w = k - l.
    // pivot is the index where the NEXT hash will be written (= current leftmost,
    // the l-mer that will be evicted on the next advance).
    //
    // Invariant when pivot = p (0 < p < w):
    //   H[p]   = leftmost l-mer of the current window (will be evicted next)
    //   H[p+1] = rightmost l-mer (just written)
    //   H[p+2..w], H[0..p-1] = the w-1 intermediate l-mers
    //
    // After reset_windows() (pivot reset to w):
    //   H[w]   = leftmost l-mer  (will be evicted on the first advance)
    //   H[0]   = rightmost l-mer (most recently written before the reset)
    //
    // M_pre[i] = min(H[0..i])   — prefix minimums from the rightmost end.
    // M_suf    = running min of the suffix (rightmost/newer) half since last reset.
    // pivot    = index into H[] where the next new hash will be written.
    //
    // hash() = min(M_pre[pivot], M_suf).
    // When pivot reaches 0 the suffix becomes the new prefix (reset_windows).

    uint64_t    H[k];
    uint64_t    M_pre[k];
    uint64_t    M_suf  = 0;
    std::size_t pivot  = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    void reset_windows() noexcept {
        const std::size_t w = k - l_;
        M_pre[0] = H[0];
        for (std::size_t i = 1; i <= w; ++i)
            M_pre[i] = umin(M_pre[i - 1], H[i]);
        M_suf  = U64_MAX;
        pivot  = w;
    }

    // Core slide: `in_2bit` is nt-encoded (A=0,C=1,T=2,G=3).
    void advance_impl(uint8_t in_2bit) noexcept {
        const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (l_ - 1))) & 0x3u;
        lmer_     = ((lmer_     & clear_msn_) << 2) | in_2bit;
        lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l_ - 1)));
        hasher_.roll(out_2bit, in_2bit);

        H[pivot] = hasher_.canonical();
        M_suf    = umin(M_suf, H[pivot]);

        if (pivot > 0)
            --pivot;
        else
            reset_windows();
    }

public:

    explicit MinimizerWindow(uint16_t l) noexcept
        : l_(l)
        , clear_msn_(~(uint64_t(0b11) << (2 * (l - 1))))
        , hasher_(l)
    {
        assert(l < k);
    }

    // Initialise to the first k characters of `seq` (ASCII ACGT, any case).
    // `seq` must contain at least k valid DNA characters.
    void reset(const char* seq) noexcept {
        lmer_     = 0;
        lmer_bar_ = 0;
        for (uint16_t i = 0; i < l_; ++i) {
            const uint8_t b = nt_hash::to_2bit(seq[i]);
            lmer_     |= (uint64_t(b)       << (2 * (l_ - 1 - i)));
            lmer_bar_ |= (uint64_t(b ^ 2u)  << (2 * i));
        }
        hasher_.init(seq);

        pivot    = k - l_;
        H[pivot] = hasher_.canonical();

        for (uint16_t i = l_; i < k; ++i) {
            const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (l_ - 1))) & 0x3u;
            const uint8_t in_2bit  = nt_hash::to_2bit(seq[i]);
            lmer_     = ((lmer_     & clear_msn_) << 2) | in_2bit;
            lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l_ - 1)));
            hasher_.roll(out_2bit, in_2bit);
            H[--pivot] = hasher_.canonical();
        }
        reset_windows();
    }

    // Slide by one ASCII character.
    void advance(char ch) noexcept {
        advance_impl(nt_hash::to_2bit(ch));
    }

    // Slide by one base in kache encoding (A=0, C=1, G=2, T=3).
    // Converts to nt encoding via b ^ (b >> 1) before processing.
    void advance_kache(uint8_t kache_b) noexcept {
        advance_impl(kache_b ^ (kache_b >> 1));
    }

    // Canonical ntHash of the l-minimizer of the current k-mer window.
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
