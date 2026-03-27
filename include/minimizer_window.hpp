#pragma once

// MinimizerWindow<k, m> — streaming m-minimizer iterator using canonical ntHash.
//
// Used in:
//   Phase 0+1 (partition_hash.hpp, partition_kmtricks.hpp) — ASCII input path.
//   Phase 2   (kache-hash Kmer_Window) — field replaces Min_Iterator<k,m> (XXH3),
//             but bypassed for hash/kmtricks superkmers (min_pos byte used instead);
//             only active for KMC-mode superkmers (min_pos == 0xFF fallback).
//
// Algorithm: two-stack (prefix/suffix) sliding window minimum over ntHash values.
// Window size w = k - m; arrays sized exactly w+1 for compact cache footprint.
// O(1) amortised per advance() call.
//
// Two input encodings are supported:
//   ASCII  — advance(char ch)         uses nt_hash::to_2bit() (A=0,C=1,T=2,G=3)
//   kache  — advance_kache(uint8_t b) accepts kache encoding (A=0,C=1,G=2,T=3),
//            converts to nt via b ^ (b >> 1), then calls advance_impl().
//
// Public interface:
//   MinimizerWindow<k, m> win;
//   win.reset(seq);              // initialise to first k ASCII characters
//   uint64_t h = win.hash();     // canonical ntHash of the m-minimizer
//   win.advance(ch);             // slide by one ASCII character
//   win.advance_kache(b);        // slide by one kache-encoded base
//   win.hash(min_coord);         // hash + dummy min_coord (for kache-hash API)
//   win.min_lmer_pos();          // absolute position of the minimizer m-mer
//                                //   (0 = first m-mer passed to reset())

#include "nt_hash.hpp"

#include <cstdint>
#include <limits>

template <uint16_t k, uint16_t m>
class MinimizerWindow {
    static_assert(m < k, "minimizer length m must be strictly less than k-mer length k");

    // Window size: number of m-mers in a k-mer window.
    static constexpr uint16_t w = k - m;

    // Mask to clear the top 2 bits of the 2*m-bit m-mer (used in advance_impl).
    static constexpr uint64_t clear_msn_ = ~(uint64_t(0b11) << (2 * (m - 1)));

    // 2-bit packed forward m-mer, used to extract the outgoing base cheaply
    // when rolling the ntHash.  nt encoding (A=0,C=1,T=2,G=3).
    // Note: lmer_bar_ (RC m-mer) was removed — roll() derives the RC contribution
    // directly from out_2bit/in_2bit via the REV[] seed table, so lmer_bar_ was
    // written every advance but never read (dead store).
    uint64_t lmer_ = 0;

    nt_hash::Roller<m> hasher_;

    // ── Two-stack sliding window minimum ─────────────────────────────────────
    //
    // H[0..w] is a circular overwrite buffer of w+1 entries.
    // pivot is the index where the NEXT hash will be written (= current leftmost,
    // the m-mer that will be evicted on the next advance).
    //
    // Invariant when pivot = p (0 < p < w):
    //   H[p]   = leftmost m-mer of the current window (will be evicted next)
    //   H[p+1] = rightmost m-mer (just written)
    //   H[p+2..w], H[0..p-1] = the w-1 intermediate m-mers
    //
    // After reset_windows() (pivot reset to w):
    //   H[w]   = leftmost m-mer  (will be evicted on the first advance)
    //   H[0]   = rightmost m-mer (most recently written before the reset)
    //
    // M_pre[i] = min(H[0..i])   — prefix minimums from the rightmost end.
    // M_suf    = running min of the suffix (rightmost/newer) half since last reset.
    // pivot    = index into H[] where the next new hash will be written.
    //
    // hash() = min(M_pre[pivot], M_suf).
    // When pivot reaches 0 the suffix becomes the new prefix (reset_windows).
    //
    // Position tracking — lazy reconstruction (not stored per base):
    //   Instead of storing H_pos[], M_pre_pos[], M_suf_pos as full uint64_t arrays
    //   in the hot loop (184 bytes written every advance), we store compact indices
    //   (uint8_t) and reconstruct absolute positions only in min_lmer_pos().
    //
    //   reset_h0_pos_   = absolute position of H[0] at the last reset_windows call.
    //   M_pre_idx[i]    = index j into H[] such that M_pre[i] == H[j]  (j ≤ i).
    //   M_suf_idx       = index j into H[] such that M_suf  == H[j]  (j > pivot).
    //
    //   next_pos_ is eliminated: instead of ++next_pos_ every base, reset_h0_pos_ is
    //   incremented by w+1 only when pivot wraps to 0 (every w+1 iterations).
    //
    //   Position reconstruction in min_lmer_pos():
    //     H[j] with j ≤ pivot (previous cycle) has position = reset_h0_pos_ - j.
    //     H[j] with j > pivot (current  cycle) has position = reset_h0_pos_ + w + 1 - j.
    //     (pivot cancels: next_pos_ + pivot - j = (reset_h0_pos_ + w - pivot + 1) + pivot - j)

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
            // Branchless select: CMOV instead of branch (50%-taken comparisons
            // cause mispredictions on every other reset_windows call otherwise).
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
        // Branchless M_suf update: avoid mispredictions on the ~21%-taken condition.
        // M_suf_idx is size_t (not uint8_t) so the compiler can emit CMOV — x86 has no
        // 8-bit CMOV, so uint8_t forced a branch here despite the ternary operator.
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
    // Separates hash rolling from window-minimum update so each phase can run
    // in a tight loop with low register pressure.  The combined advance_kache()
    // keeps fwd_/rev_ live alongside all the partition-logic variables (pivot,
    // M_suf, sk_start, prev_hash, pid, …), causing the compiler to spill them
    // to the stack and incur store-to-load forwarding latency on the loop-carried
    // dependency chain.
    //
    // Two-pass usage (equivalent to N consecutive advance_kache calls):
    //   for (size_t i = 0; i < N; i++) hash_buf[i] = win.roll_hash_kache(seq[i]);
    //   for (size_t i = 0; i < N; i++) win.advance_with_hash(hash_buf[i]);
    //
    // Pass 1 only keeps {fwd_, rev_, lmer_, ptr, counter} live — all fit in GPRs,
    // eliminating the SLF bottleneck on fwd_.
    // Pass 2 processes the pre-computed hash buffer without touching the hasher.

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

    // Absolute position (0-indexed from the sequence passed to reset()) of the
    // m-mer that achieves the minimizer hash in the current k-mer window.
    // Consistent with hash(): same tie-breaking (suffix wins on equality).
    //
    // Positions are reconstructed lazily (not stored per base):
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
