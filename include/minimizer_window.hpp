#pragma once

// MinimizerWindow<k> — streaming l-minimizer iterator using canonical ntHash.
//
// Drop-in replacement for cuttlefish::Min_Iterator<k> in Phase 1 (partition_hash.hpp)
// and Phase 0 (partition_kmtricks.hpp).
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ntHash values.
// O(k) space, O(1) amortised per advance() call, identical to the data structure
// used by cuttlefish::Min_Iterator<k> — only the hash function changes:
//   OLD: min(wyhash(lmer), wyhash(lmer_bar))
//   NEW: nt_hash::Roller::canonical() = fwd_hash XOR rev_hash
//        (same convention as simd-minimizers, seq-hash crate, NtHasher<true>)
//
// Public interface:
//   MinimizerWindow<k> win(l);
//   win.reset(seq);           // initialise to first k characters
//   uint64_t h = win.hash();  // canonical ntHash of the l-minimizer
//   win.advance(ch);          // slide by one character
//   h = win.hash();

#include "nt_hash.hpp"

#include <cstdint>
#include <cassert>
#include <limits>

template <uint16_t k>
class MinimizerWindow {

    uint16_t l_;
    const uint64_t clear_msn_;  // mask to clear the top 2 bits of the 2*l-bit l-mer

    // 2-bit packed l-mer (forward and RC), tracked to extract the outgoing
    // base cheaply when rolling the ntHash.
    uint64_t lmer_     = 0;
    uint64_t lmer_bar_ = 0;

    nt_hash::Roller hasher_;

    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w] stores canonical l-mer hashes in REVERSED position order:
    //   H[w-1] = hash of the FIRST (leftmost) l-mer in the current k-mer
    //   H[0]   = hash of the LAST  (rightmost) l-mer in the current k-mer
    // where w = k - l (= number of l-mer positions per k-mer minus 1).
    //
    // M_pre[i] = min(H[0..i])   — prefix minimums from the right end.
    // M_suf    = running min of the current suffix (left) half.
    // pivot    = index into H[] where the next new hash will be written.
    //
    // hash() = min(M_pre[pivot], M_suf).
    // When pivot wraps to 0 the suffix becomes the new prefix (reset_windows).

    uint64_t    H[k];
    uint64_t    M_pre[k];
    uint64_t    M_suf  = 0;
    std::size_t pivot  = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();
    static constexpr auto umin = [](uint64_t a, uint64_t b) noexcept { return a < b ? a : b; };

    // Recompute prefix mins from the current H array, reset M_suf, restore pivot.
    // Called when the suffix half has been filled (pivot reached 0).
    void reset_windows() noexcept {
        const std::size_t w = k - l_;
        M_pre[0] = H[0];
        for (std::size_t i = 1; i <= w; ++i)
            M_pre[i] = umin(M_pre[i - 1], H[i]);
        M_suf  = U64_MAX;
        pivot  = w;
    }

public:

    explicit MinimizerWindow(uint16_t l) noexcept
        : l_(l)
        , clear_msn_(~(uint64_t(0b11) << (2 * (l - 1))))
        , hasher_(l)
    {
        assert(l < k);
    }

    // Initialise to the first k characters of `seq`.
    // `seq` must contain at least k ACGT characters (no N, no newlines).
    void reset(const char* seq) noexcept {
        // Build the first l-mer from seq[0..l-1].
        lmer_     = 0;
        lmer_bar_ = 0;
        for (uint16_t i = 0; i < l_; ++i) {
            const uint8_t b = nt_hash::to_2bit(seq[i]);
            lmer_     |= (uint64_t(b)               << (2 * (l_ - 1 - i)));
            lmer_bar_ |= (uint64_t(b ^ 2u)           << (2 * i));
        }
        hasher_.init(seq);

        // Fill H in reversed order: H[w] = hash(lmer@0), …, H[0] = hash(lmer@w).
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
        // pivot == 0 after the loop; set up prefix-min arrays.
        reset_windows();
    }

    // Slide the window one position to the right by consuming `ch`.
    // `ch` must be an ACGT character.
    void advance(char ch) noexcept {
        const uint8_t out_2bit = static_cast<uint8_t>(lmer_ >> (2 * (l_ - 1))) & 0x3u;
        const uint8_t in_2bit  = nt_hash::to_2bit(ch);
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

    // Canonical ntHash of the l-minimizer of the current k-mer window.
    uint64_t hash() const noexcept {
        return umin(M_pre[pivot], M_suf);
    }
};
