#pragma once

// Portable 256-bit SIMD abstraction used by kache-hash.
// Provides the five operations needed to vectorise 32-byte checksum arrays:
//   simd256_t              — register type
//   simd256_load(ptr)      — aligned/unaligned 32-byte load
//   simd256_zero()         — all-zeros register
//   simd256_set1(c)        — broadcast one byte to all 32 positions
//   simd256_cmpeq(a, b)    — per-byte equality (0xFF match / 0x00 no-match)
//   simd256_movemask(v)    — collapse MSB of each byte → uint32_t bitmask
//
// x86-64  : AVX2 (single 256-bit register, one instruction per op)
// AArch64 : NEON (two 128-bit registers paired in a struct)

#include <cstdint>

// ── x86-64 / AVX2 ────────────────────────────────────────────────────────────
#if defined(__x86_64__) || defined(_M_X64)

#include <immintrin.h>

using simd256_t = __m256i;

inline simd256_t simd256_load(const void* p) noexcept
    { return _mm256_load_si256(static_cast<const __m256i*>(p)); }

inline simd256_t simd256_zero() noexcept
    { return _mm256_setzero_si256(); }

inline simd256_t simd256_set1(char c) noexcept
    { return _mm256_set1_epi8(c); }

inline simd256_t simd256_cmpeq(simd256_t a, simd256_t b) noexcept
    { return _mm256_cmpeq_epi8(a, b); }

inline uint32_t simd256_movemask(simd256_t v) noexcept
    { return static_cast<uint32_t>(_mm256_movemask_epi8(v)); }


// ── AArch64 / NEON ───────────────────────────────────────────────────────────
#elif defined(__aarch64__) || defined(_M_ARM64)

#include <arm_neon.h>

struct simd256_t { uint8x16_t lo, hi; };

inline simd256_t simd256_load(const void* p) noexcept {
    const auto* b = static_cast<const uint8_t*>(p);
    return { vld1q_u8(b), vld1q_u8(b + 16) };
}

inline simd256_t simd256_zero() noexcept
    { return { vdupq_n_u8(0), vdupq_n_u8(0) }; }

inline simd256_t simd256_set1(char c) noexcept {
    const auto v = vdupq_n_u8(static_cast<uint8_t>(c));
    return { v, v };
}

inline simd256_t simd256_cmpeq(simd256_t a, simd256_t b) noexcept
    { return { vceqq_u8(a.lo, b.lo), vceqq_u8(a.hi, b.hi) }; }

// Collapse the MSB of each byte in a 16-byte NEON register into a 16-bit mask.
// Classic sse2neon technique: shift each MSB to its bit-index position, then
// three rounds of pairwise-add to collect all 8 bits into one byte.
static inline uint16_t neon_movemask16(uint8x16_t v) noexcept {
    static const int8_t kShifts[8] = { -7, -6, -5, -4, -3, -2, -1, 0 };
    const int8x8_t sv = vld1_s8(kShifts);
    uint8x8_t lo = vand_u8(vget_low_u8(v),  vdup_n_u8(0x80));
    uint8x8_t hi = vand_u8(vget_high_u8(v), vdup_n_u8(0x80));
    lo = vreinterpret_u8_s8(vshl_s8(vreinterpret_s8_u8(lo), sv));
    hi = vreinterpret_u8_s8(vshl_s8(vreinterpret_s8_u8(hi), sv));
    uint8x8_t r = vpadd_u8(lo, hi);
    r = vpadd_u8(r, r);
    r = vpadd_u8(r, r);
    return static_cast<uint16_t>(vget_lane_u8(r, 0))
         | (static_cast<uint16_t>(vget_lane_u8(r, 1)) << 8);
}

inline uint32_t simd256_movemask(simd256_t v) noexcept {
    return static_cast<uint32_t>(neon_movemask16(v.lo))
         | (static_cast<uint32_t>(neon_movemask16(v.hi)) << 16);
}


#else
#error "Unsupported architecture: tuna requires x86-64 (AVX2) or AArch64 (NEON)"
#endif
