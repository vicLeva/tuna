#pragma once

// MinimizerWindow<k, m> — streaming m-minimizer iterator using canonical ntHash.
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ntHash values.
// Window size w = k - m; arrays sized exactly w+1. O(1) amortised per advance().
//
// Two input encodings:
//   ASCII  — advance(char ch)         (to_2bit: A=0,C=1,T=2,G=3)
//   kache  — advance_kache(uint8_t b) (A=0,C=1,G=2,T=3, converted via b^(b>>1))
//
// Usage:
//   MinimizerWindow<k, m> win;
//   win.reset(seq);          // init from first k ASCII characters
//   win.advance(ch);         // slide by one character
//   win.hash();              // canonical ntHash of the current window's minimizer
//   win.min_lmer_pos();      // absolute position of the minimizer (0 = first m-mer)

#include "nt_hash.hpp"

#include <cstdint>
#include <limits>

template <uint16_t k, uint16_t m>
class MinimizerWindow {
    static_assert(m < k,  "minimizer length m must be strictly less than k-mer length k");
    // lmer_ stores the 2-bit packed m-mer in a uint64_t — limits m to 32.
    static_assert(m <= 32, "minimizer length m must be <= 32 (lmer_ is uint64_t; 2*m bits required)");

    // Window size: number of m-mers in a k-mer window.
    static constexpr uint16_t w = k - m;

    // Mask to clear the top 2 bits of the 2*m-bit m-mer (used in advance_impl).
    static constexpr uint64_t clear_msn_ = ~(uint64_t(0b11) << (2 * (m - 1)));

    // 2-bit packed forward m-mer in nt encoding (A=0,C=1,T=2,G=3).
    // Used to extract the outgoing base cheaply when rolling the ntHash.
    uint64_t lmer_ = 0;

    nt_hash::Roller<m> hasher_;

    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w]: circular buffer of w+1 hash values. pivot = next write index
    //   (= leftmost m-mer, evicted on the next advance).
    // M_pre[i] = min(H[0..i])  — prefix minimums from the rightmost end.
    // M_suf    = running min of hashes written since the last reset_windows().
    // hash() = min(M_pre[pivot], M_suf).  pivot reaches 0 → reset_windows().
    //
    // Position tracking: M_pre_idx[i] / M_suf_idx store indices into H[] rather
    // than full uint64_t positions (saves 184 bytes of hot-loop writes).
    // reset_h0_pos_ = absolute position of H[0] at the last reset_windows();
    // positions are reconstructed lazily in min_lmer_pos().

    uint64_t    H[w + 1];
    uint64_t    M_pre[w + 1];
    uint8_t     M_pre_idx[w + 1] = {}; // index into H[] for each prefix-min entry
    uint64_t    M_suf      = 0;
    std::size_t M_suf_idx  = 0;        // index into H[] for the current suffix min
    uint64_t    reset_h0_pos_ = 0;     // position of H[0] at the last reset_windows()
    std::size_t pivot      = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    void reset_windows() noexcept {
        M_pre[0]     = H[0];
        M_pre_idx[0] = 0;
        for (std::size_t i = 1; i <= w; ++i) {
            // Branchless CMOV: 50%-taken comparisons cause mispredictions otherwise.
            const bool lt = H[i] < M_pre[i - 1];
            M_pre[i]     = lt ? H[i]                    : M_pre[i - 1];
            M_pre_idx[i] = lt ? static_cast<uint8_t>(i) : M_pre_idx[i - 1];
        }
        M_suf     = U64_MAX;
        M_suf_idx = 0;
        pivot     = w;
    }

    // Core slide: `in_2bit` is nt-encoded (A=0,C=1,T=2,G=3).
    void advance_impl(uint8_t in_2bit) noexcept {
        const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (m - 1))) & 0x3u;
        lmer_ = ((lmer_ & clear_msn_) << 2) | in_2bit;
        hasher_.roll(out_2bit, in_2bit);

        const uint64_t h = hasher_.canonical();
        H[pivot] = h;
        // Branchless CMOV: M_suf_idx is size_t (not uint8_t) — x86 has no 8-bit CMOV.
        const bool lt = (h < M_suf);
        M_suf     = lt ? h      : M_suf;
        M_suf_idx = lt ? pivot  : M_suf_idx;

        if (pivot > 0)
            --pivot;
        else {
            reset_h0_pos_ += static_cast<uint64_t>(w) + 1; // advance by one full cycle
            reset_windows();
        }
    }

public:

    MinimizerWindow() noexcept = default;

    // Initialise to the first k characters of `seq` (ASCII ACGT, any case).
    // `seq` must contain at least k valid DNA characters.
    // After reset(), positions are 0-indexed from seq[0].
    void reset(const char* seq) noexcept {
        lmer_ = 0;
        for (uint16_t i = 0; i < m; ++i)
            lmer_ |= (uint64_t(nt_hash::to_2bit(seq[i])) << (2 * (m - 1 - i)));
        hasher_.init(seq);

        pivot    = w;
        H[pivot] = hasher_.canonical(); // m-mer seq[0..m-1] at position 0

        for (uint16_t i = m; i < k; ++i) {
            const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (m - 1))) & 0x3u;
            const uint8_t in_2bit  = nt_hash::to_2bit(seq[i]);
            lmer_ = ((lmer_ & clear_msn_) << 2) | in_2bit;
            hasher_.roll(out_2bit, in_2bit);
            --pivot;
            H[pivot] = hasher_.canonical(); // m-mer seq[i-m+1..i] at position i-m+1
        }
        // After loop: pivot=0, H[0] = hash of last m-mer in the initial k-mer window
        // (at position w = k-m).  First advance gives position w+1.
        reset_h0_pos_ = static_cast<uint64_t>(w); // position of H[0] in this initial fill
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

    // ── Two-phase advance API ─────────────────────────────────────────────────
    //
    // Separates hash rolling (roll_hash_kache) from window-min update
    // (advance_with_hash) to allow each phase to run with lower register pressure.
    // Equivalent to N consecutive advance_kache() calls when used in tandem.

    // Roll the hasher by one kache-encoded base; return the canonical hash.
    // Precondition: advance_with_hash() must be called with the returned value
    //               before any hash()/min_lmer_pos() query.
    uint64_t roll_hash_kache(uint8_t kache_b) noexcept {
        const uint8_t in_2bit  = kache_b ^ (kache_b >> 1);
        const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (m - 1))) & 0x3u;
        lmer_ = ((lmer_ & clear_msn_) << 2) | in_2bit;
        hasher_.roll(out_2bit, in_2bit);
        return hasher_.canonical();
    }

    // Advance the sliding-window minimum with a pre-computed hash value.
    // Must be called with the value returned by roll_hash_kache(), same order.
    void advance_with_hash(uint64_t h) noexcept {
        H[pivot] = h;
        const bool lt = (h < M_suf);
        M_suf     = lt ? h      : M_suf;
        M_suf_idx = lt ? pivot  : M_suf_idx;

        if (pivot > 0)
            --pivot;
        else {
            reset_h0_pos_ += static_cast<uint64_t>(w) + 1;
            reset_windows();
        }
    }

    // Canonical ntHash of the m-minimizer of the current k-mer window.
    uint64_t hash() const noexcept {
        return umin(M_pre[pivot], M_suf);
    }

    // Same as hash(), but also sets min_coord = 0 (position/orientation not
    // tracked; provided for compatibility with the kache-hash Kmer_Window API).
    uint64_t hash(uint8_t& min_coord) const noexcept {
        min_coord = 0;
        return hash();
    }

    // Absolute position of the minimizer m-mer (0 = first m-mer passed to reset()).
    // Reconstructed lazily from stored H[] indices:
    //   j ≤ pivot (previous cycle): position = reset_h0_pos_ - j
    //   j > pivot (current  cycle): position = reset_h0_pos_ + w + 1 - j
    uint64_t min_lmer_pos() const noexcept {
        if (M_pre[pivot] < M_suf) {
            const std::size_t j = M_pre_idx[pivot]; // j ≤ pivot (previous cycle)
            return reset_h0_pos_ - j;
        } else {
            const std::size_t j = M_suf_idx;        // j > pivot (current cycle)
            return reset_h0_pos_ + static_cast<uint64_t>(w) + 1 - j;
        }
    }
};
