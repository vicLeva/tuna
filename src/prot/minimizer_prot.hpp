#pragma once

// MinimizerWindowProt<k, m> — forward-only streaming m-minimizer for proteins.
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ProtRoller
// hash values. Window size w = k - m; O(1) amortised per advance().
//
// lmer_ stores the last m amino acid codes (5 bits each) in a uint64_t:
//   bits[5*(m-1)+4 : 5*(m-1)] = oldest (outgoing on next advance)
//   bits[4:0]                  = newest
// Valid for m ≤ 12 (5*12 = 60 bits ≤ 64).
//
// Usage:
//   MinimizerWindowProt<k, m> win;
//   win.reset(enc);        // init from first k encoded amino acids
//   win.advance(in5);      // slide by one encoded amino acid
//   win.hash();            // hash of current window's minimizer m-mer
//   win.min_lmer_pos();    // absolute position of the minimizer

#include "aa_hash.hpp"
#include <cstdint>
#include <limits>

namespace prot {

template <uint16_t k, uint16_t m>
class MinimizerWindowProt {
    static_assert(m < k,   "minimizer m must be strictly less than k");
    static_assert(m >= 2,  "minimizer m must be >= 2");
    static_assert(m <= 12, "minimizer m must be <= 12 (5*m bits must fit in uint64_t)");

    static constexpr uint16_t w = k - m;  // window = number of distinct m-mers in a k-mer

    // Mask for the 5*m least-significant bits of lmer_.
    static constexpr uint64_t lmer_mask =
        (m == 12) ? uint64_t(0x0FFFFFFFFFFFFFFF)
                  : ((uint64_t(1) << (5 * m)) - 1);
    // Shift to extract the outgoing (oldest) amino acid from lmer_.
    static constexpr int outgoing_shift = 5 * (m - 1);

    uint64_t      lmer_ = 0;
    ProtRoller<m> hasher_;

    // Two-stack sliding window minimum buffers.
    uint64_t    H[w + 1];
    uint64_t    M_pre[w + 1];
    uint8_t     M_pre_idx[w + 1] = {};
    uint64_t    M_suf      = 0;
    size_t      M_suf_idx  = 0;
    uint64_t    reset_h0_pos_ = 0;
    size_t      pivot = 0;

    static constexpr uint64_t U64_MAX = std::numeric_limits<uint64_t>::max();

    void reset_windows() noexcept {
        M_pre[0]     = H[0];
        M_pre_idx[0] = 0;
        for (size_t i = 1; i <= w; ++i) {
            const bool lt    = H[i] < M_pre[i - 1];
            M_pre[i]         = lt ? H[i]                    : M_pre[i - 1];
            M_pre_idx[i]     = lt ? static_cast<uint8_t>(i) : M_pre_idx[i - 1];
        }
        M_suf     = U64_MAX;
        M_suf_idx = 0;
        pivot     = w;
    }

    void advance_impl(uint8_t in5) noexcept {
        const uint8_t out5 = static_cast<uint8_t>(lmer_ >> outgoing_shift) & 0x1fu;
        lmer_ = ((lmer_ << 5) | in5) & lmer_mask;
        hasher_.roll(out5, in5);
        const uint64_t h = hasher_.hash();
        H[pivot] = h;
        const bool lt = (h < M_suf);
        M_suf     = lt ? h     : M_suf;
        M_suf_idx = lt ? pivot : M_suf_idx;
        if (pivot > 0)
            --pivot;
        else {
            reset_h0_pos_ += static_cast<uint64_t>(w) + 1;
            reset_windows();
        }
    }

public:
    MinimizerWindowProt() noexcept = default;

    // Initialise from the first k encoded amino acids (5-bit values, 0-19).
    void reset(const uint8_t* enc) noexcept {
        // Pack enc[0..m-1] into lmer_, oldest at top.
        lmer_ = 0;
        for (uint16_t i = 0; i < m; ++i)
            lmer_ = ((lmer_ << 5) | enc[i]) & lmer_mask;

        hasher_.init(enc);
        pivot    = w;
        H[pivot] = hasher_.hash();   // m-mer at position 0

        for (uint16_t i = m; i < k; ++i) {
            const uint8_t out5 = static_cast<uint8_t>(lmer_ >> outgoing_shift) & 0x1fu;
            const uint8_t in5  = enc[i];
            lmer_ = ((lmer_ << 5) | in5) & lmer_mask;
            hasher_.roll(out5, in5);
            --pivot;
            H[pivot] = hasher_.hash();
        }
        reset_h0_pos_ = static_cast<uint64_t>(w);
        reset_windows();
    }

    // Slide by one encoded amino acid (5-bit value 0-19).
    void advance(uint8_t in5) noexcept { advance_impl(in5); }

    // Hash of the minimizer m-mer of the current k-mer window.
    uint64_t hash() const noexcept {
        return M_pre[pivot] < M_suf ? M_pre[pivot] : M_suf;
    }

    // Same as hash(), writes min_coord=0 (provided for API compatibility).
    uint64_t hash(uint8_t& min_coord) const noexcept {
        min_coord = 0;
        return hash();
    }

    // Absolute position of the minimizer m-mer (0 = first m-mer passed to reset()).
    uint64_t min_lmer_pos() const noexcept {
        if (M_pre[pivot] < M_suf) {
            const size_t j = M_pre_idx[pivot];
            return reset_h0_pos_ - j;
        } else {
            const size_t j = M_suf_idx;
            return reset_h0_pos_ + static_cast<uint64_t>(w) + 1 - j;
        }
    }
};

} // namespace prot
