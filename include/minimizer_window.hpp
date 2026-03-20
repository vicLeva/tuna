#pragma once

// MinimizerWindow<k, l>    — streaming l-minimizer iterator using canonical ntHash.
// MinimizerWindowXXH<k, l> — same algorithm, but uses canonical XXH3_64bits on the
//                             2-bit packed l-mer instead of ntHash rolling.
//                             Used by the -xxhash partition strategy benchmark.
//
// Used in:
//   Phase 0+1 (partition_hash.hpp, partition_kmtricks.hpp) — ASCII input path.
//   Phase 2   (kache-hash Kmer_Window) — field replaces Min_Iterator<k,l> (XXH3),
//             but bypassed for hash/kmtricks superkmers (min_pos byte used instead);
//             only active for KMC-mode superkmers (min_pos == 0xFF fallback).
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ntHash values.
// Window size w = k - l; arrays sized exactly w+1 for compact cache footprint.
// O(1) amortised per advance() call.
//
// Two input encodings are supported:
//   ASCII  — advance(char ch)         uses nt_hash::to_2bit() (A=0,C=1,T=2,G=3)
//   kache  — advance_kache(uint8_t b) accepts kache encoding (A=0,C=1,G=2,T=3),
//            converts to nt via b ^ (b >> 1), then calls advance_impl().
//
// Public interface:
//   MinimizerWindow<k, l> win;
//   win.reset(seq);              // initialise to first k ASCII characters
//   uint64_t h = win.hash();     // canonical ntHash of the l-minimizer
//   win.advance(ch);             // slide by one ASCII character
//   win.advance_kache(b);        // slide by one kache-encoded base
//   win.hash(min_coord);         // hash + dummy min_coord (for kache-hash API)
//   win.min_lmer_pos();          // absolute position of the minimizer l-mer
//                                //   (0 = first l-mer passed to reset())

#include "nt_hash.hpp"
#include "xxHash/xxhash.h"

#include <cstdint>
#include <limits>

template <uint16_t k, uint16_t l>
class MinimizerWindow {
    static_assert(l < k, "minimizer length l must be strictly less than k-mer length k");

    // Window size: number of l-mers in a k-mer window.
    static constexpr uint16_t w = k - l;

    // Mask to clear the top 2 bits of the 2*l-bit l-mer (used in advance_impl).
    static constexpr uint64_t clear_msn_ = ~(uint64_t(0b11) << (2 * (l - 1)));

    // 2-bit packed l-mer (forward and RC), tracked to extract the outgoing
    // base cheaply when rolling the ntHash.  nt encoding (A=0,C=1,T=2,G=3).
    uint64_t lmer_     = 0;
    uint64_t lmer_bar_ = 0;

    nt_hash::Roller<l> hasher_;

    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w] is a circular overwrite buffer of w+1 entries.
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
    //
    // Position tracking (H_pos[], M_pre_pos[], M_suf_pos):
    //   Each H[i] has a paired H_pos[i] = absolute position (0-indexed from
    //   the first l-mer passed to reset()) of that l-mer.  M_pre_pos[i] and
    //   M_suf_pos mirror the position of the minimum alongside the hash.
    //   min_lmer_pos() returns the position of the current minimizer with no
    //   extra work — avoids the O(superkmer_len) find_minimizer_pos rescan.

    uint64_t    H[w + 1];
    uint64_t    H_pos[w + 1];     // absolute l-mer position for each H[] entry
    uint64_t    M_pre[w + 1];
    uint64_t    M_pre_pos[w + 1]; // position of the minimum in each prefix slice
    uint64_t    M_suf     = 0;
    uint64_t    M_suf_pos = 0;
    std::size_t pivot     = 0;
    uint64_t    next_pos_ = 0;    // absolute position of the next l-mer to be written

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    void reset_windows() noexcept {
        M_pre[0]     = H[0];
        M_pre_pos[0] = H_pos[0];
        for (std::size_t i = 1; i <= w; ++i) {
            if (H[i] < M_pre[i - 1]) {
                M_pre[i]     = H[i];
                M_pre_pos[i] = H_pos[i];
            } else {
                M_pre[i]     = M_pre[i - 1];
                M_pre_pos[i] = M_pre_pos[i - 1];
            }
        }
        M_suf     = U64_MAX;
        M_suf_pos = 0;
        pivot     = w;
    }

    // Core slide: `in_2bit` is nt-encoded (A=0,C=1,T=2,G=3).
    void advance_impl(uint8_t in_2bit) noexcept {
        const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (l - 1))) & 0x3u;
        lmer_     = ((lmer_     & clear_msn_) << 2) | in_2bit;
        lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l - 1)));
        hasher_.roll(out_2bit, in_2bit);

        const uint64_t h    = hasher_.canonical();
        const uint64_t pos  = next_pos_++;
        H[pivot]     = h;
        H_pos[pivot] = pos;
        if (h < M_suf) { M_suf = h; M_suf_pos = pos; }

        if (pivot > 0)
            --pivot;
        else
            reset_windows();
    }

public:

    MinimizerWindow() noexcept = default;

    // Initialise to the first k characters of `seq` (ASCII ACGT, any case).
    // `seq` must contain at least k valid DNA characters.
    // After reset(), positions are 0-indexed from seq[0].
    void reset(const char* seq) noexcept {
        lmer_     = 0;
        lmer_bar_ = 0;
        for (uint16_t i = 0; i < l; ++i) {
            const uint8_t b = nt_hash::to_2bit(seq[i]);
            lmer_     |= (uint64_t(b)      << (2 * (l - 1 - i)));
            lmer_bar_ |= (uint64_t(b ^ 2u) << (2 * i));
        }
        hasher_.init(seq);

        pivot       = w;
        H[pivot]    = hasher_.canonical();
        H_pos[pivot]= 0; // l-mer seq[0..l-1] is at position 0

        for (uint16_t i = l; i < k; ++i) {
            const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (l - 1))) & 0x3u;
            const uint8_t in_2bit  = nt_hash::to_2bit(seq[i]);
            lmer_     = ((lmer_     & clear_msn_) << 2) | in_2bit;
            lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l - 1)));
            hasher_.roll(out_2bit, in_2bit);
            --pivot;
            H[pivot]     = hasher_.canonical();
            H_pos[pivot] = static_cast<uint64_t>(i - l + 1); // l-mer seq[i-l+1..i]
        }
        next_pos_ = static_cast<uint64_t>(w + 1); // next advance adds l-mer at position w+1
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

    // Absolute position (0-indexed from the sequence passed to reset()) of the
    // l-mer that achieves the minimizer hash in the current k-mer window.
    // Consistent with hash(): same tie-breaking (suffix wins on equality).
    uint64_t min_lmer_pos() const noexcept {
        return (M_pre[pivot] < M_suf) ? M_pre_pos[pivot] : M_suf_pos;
    }
};


// ─── MinimizerWindowXXH<k, l> ────────────────────────────────────────────────
//
// Drop-in replacement for MinimizerWindow<k, l> that hashes each l-mer with
// canonical XXH3_64bits instead of canonical ntHash.
//
// "Canonical" here mirrors the ntHash convention:
//   hash = min(XXH3_64bits(&lmer_fwd, 8), XXH3_64bits(&lmer_rc, 8))
// where lmer_fwd and lmer_rc are the 2-bit packed l-mer and its reverse
// complement, maintained identically to MinimizerWindow.
//
// The sliding-window minimum structure (two-stack prefix/suffix) is unchanged.
// The public interface is identical to MinimizerWindow.

template <uint16_t k, uint16_t l>
class MinimizerWindowXXH {
    static_assert(l < k, "minimizer length l must be strictly less than k-mer length k");

    static constexpr uint16_t  w         = k - l;
    static constexpr uint64_t  clear_msn = ~(uint64_t(0b11) << (2 * (l - 1)));

    uint64_t lmer_     = 0;   // 2-bit packed forward l-mer (nt encoding)
    uint64_t lmer_bar_ = 0;   // 2-bit packed reverse-complement l-mer

    // Two-stack sliding window minimum (same layout as MinimizerWindow).
    uint64_t    H[w + 1];
    uint64_t    H_pos[w + 1];
    uint64_t    M_pre[w + 1];
    uint64_t    M_pre_pos[w + 1];
    uint64_t    M_suf     = 0;
    uint64_t    M_suf_pos = 0;
    std::size_t pivot     = 0;
    uint64_t    next_pos_ = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    static uint64_t xxh_canonical(uint64_t fwd, uint64_t rc) noexcept {
        const uint64_t hf = XXH3_64bits(&fwd, sizeof(uint64_t));
        const uint64_t hr = XXH3_64bits(&rc,  sizeof(uint64_t));
        return hf < hr ? hf : hr;
    }

    void reset_windows() noexcept {
        M_pre[0]     = H[0];
        M_pre_pos[0] = H_pos[0];
        for (std::size_t i = 1; i <= w; ++i) {
            if (H[i] < M_pre[i - 1]) {
                M_pre[i]     = H[i];
                M_pre_pos[i] = H_pos[i];
            } else {
                M_pre[i]     = M_pre[i - 1];
                M_pre_pos[i] = M_pre_pos[i - 1];
            }
        }
        M_suf     = U64_MAX;
        M_suf_pos = 0;
        pivot     = w;
    }

    void advance_impl(uint8_t in_2bit) noexcept {
        lmer_     = ((lmer_     & clear_msn) << 2) | in_2bit;
        lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l - 1)));

        const uint64_t h   = xxh_canonical(lmer_, lmer_bar_);
        const uint64_t pos = next_pos_++;
        H[pivot]     = h;
        H_pos[pivot] = pos;
        if (h < M_suf) { M_suf = h; M_suf_pos = pos; }

        if (pivot > 0) --pivot;
        else           reset_windows();
    }

public:

    MinimizerWindowXXH() noexcept = default;

    void reset(const char* seq) noexcept {
        lmer_     = 0;
        lmer_bar_ = 0;
        for (uint16_t i = 0; i < l; ++i) {
            const uint8_t b = nt_hash::to_2bit(seq[i]);
            lmer_     |= (uint64_t(b)      << (2 * (l - 1 - i)));
            lmer_bar_ |= (uint64_t(b ^ 2u) << (2 * i));
        }

        pivot        = w;
        H[pivot]     = xxh_canonical(lmer_, lmer_bar_);
        H_pos[pivot] = 0;

        for (uint16_t i = l; i < k; ++i) {
            const uint8_t in_2bit = nt_hash::to_2bit(seq[i]);
            lmer_     = ((lmer_     & clear_msn) << 2) | in_2bit;
            lmer_bar_ = (lmer_bar_ >> 2) | (uint64_t(in_2bit ^ 2u) << (2 * (l - 1)));
            --pivot;
            H[pivot]     = xxh_canonical(lmer_, lmer_bar_);
            H_pos[pivot] = static_cast<uint64_t>(i - l + 1);
        }
        next_pos_ = static_cast<uint64_t>(w + 1);
        reset_windows();
    }

    void advance(char ch) noexcept {
        advance_impl(nt_hash::to_2bit(ch));
    }

    void advance_kache(uint8_t kache_b) noexcept {
        advance_impl(kache_b ^ (kache_b >> 1));
    }

    uint64_t hash() const noexcept {
        return umin(M_pre[pivot], M_suf);
    }

    uint64_t hash(uint8_t& min_coord) const noexcept {
        min_coord = 0;
        return hash();
    }

    uint64_t min_lmer_pos() const noexcept {
        return (M_pre[pivot] < M_suf) ? M_pre_pos[pivot] : M_suf_pos;
    }
};
