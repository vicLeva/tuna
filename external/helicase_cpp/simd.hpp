#pragma once
#include "config.hpp"
#include <cstdint>
#include <cstring>

// Forward declarations of bitmask structs (defined in lexer.hpp, but needed here).
// We define them here to avoid circular includes.
namespace helicase {

struct FastaBitmask {
    uint64_t open_bracket;
    uint64_t line_feeds;
    uint64_t is_dna;
    __uint128_t two_bits;
    uint64_t high_bit;
    uint64_t low_bit;
    uint64_t mask_non_actg;
    uint64_t mask_n;
};

struct FastqBitmask {
    uint64_t line_feeds;
    uint64_t is_dna;
    __uint128_t two_bits;
    uint64_t high_bit;
    uint64_t low_bit;
    uint64_t mask_non_actg;
    uint64_t mask_n;
};

} // namespace helicase

// ─── AVX2 path ───────────────────────────────────────────────────────────────
#if defined(__AVX2__)
#include <immintrin.h>

namespace helicase::simd {

namespace {

// LUT for ACTG detection via VPSHUFB: "A_C_T_G_________" repeated twice.
inline __m256i make_lut_actg() {
    alignas(32) static const char s[33] = "A_C_T_G_________A_C_T_G_________";
    return _mm256_load_si256((const __m256i*)s);
}

inline uint64_t movemask_64(__m256i v1, __m256i v2) {
    return ((uint64_t)(uint32_t)_mm256_movemask_epi8(v1)) |
           (((uint64_t)(uint32_t)_mm256_movemask_epi8(v2)) << 32);
}

inline uint64_t u8_mask_avx2(__m256i v1, __m256i v2, __m256i vc) {
    return movemask_64(_mm256_cmpeq_epi8(v1, vc), _mm256_cmpeq_epi8(v2, vc));
}

template <uint64_t CONFIG>
inline void bitpack_dna(__m256i v1, __m256i v2,
                        uint64_t& is_dna, __uint128_t& two_bits,
                        uint64_t& high_bit, uint64_t& low_bit,
                        uint64_t& mask_non_actg, uint64_t& mask_n)
{
    using namespace advanced;

    if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR)) {
        uint64_t mm_hi_1 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v1, 5));
        uint64_t mm_lo_1 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v1, 6));
        uint64_t mm_hi_2 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v2, 5));
        uint64_t mm_lo_2 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v2, 6));
        high_bit = mm_hi_1 | (mm_hi_2 << 32);
        low_bit  = mm_lo_1 | (mm_lo_2 << 32);
    }

    if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_PACKED)) {
#if defined(__BMI2__)
        uint64_t mm_hi_1 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v1, 5));
        uint64_t mm_lo_1 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v1, 6));
        uint64_t mm_hi_2 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v2, 5));
        uint64_t mm_lo_2 = (uint32_t)_mm256_movemask_epi8(_mm256_slli_epi16(v2, 6));
        uint64_t mm_1 = _pdep_u64(mm_hi_1, 0xAAAAAAAAAAAAAAAAull) |
                        _pdep_u64(mm_lo_1, 0x5555555555555555ull);
        uint64_t mm_2 = _pdep_u64(mm_hi_2, 0xAAAAAAAAAAAAAAAAull) |
                        _pdep_u64(mm_lo_2, 0x5555555555555555ull);
        two_bits = (__uint128_t)mm_1 | ((__uint128_t)mm_2 << 64);
#else
        // No-PDEP path: permute + unpack
        __m256i iv1 = _mm256_permute4x64_epi64(v1, 0xD8);
        __m256i iv2 = _mm256_permute4x64_epi64(v2, 0xD8);
        __m256i hi_1 = _mm256_slli_epi16(iv1, 5), lo_1 = _mm256_slli_epi16(iv1, 6);
        __m256i hi_2 = _mm256_slli_epi16(iv2, 5), lo_2 = _mm256_slli_epi16(iv2, 6);
        uint64_t mm_hi_1 = (uint32_t)_mm256_movemask_epi8(_mm256_unpackhi_epi8(lo_1, hi_1));
        uint64_t mm_lo_1 = (uint32_t)_mm256_movemask_epi8(_mm256_unpacklo_epi8(lo_1, hi_1));
        uint64_t mm_hi_2 = (uint32_t)_mm256_movemask_epi8(_mm256_unpackhi_epi8(lo_2, hi_2));
        uint64_t mm_lo_2 = (uint32_t)_mm256_movemask_epi8(_mm256_unpacklo_epi8(lo_2, hi_2));
        uint64_t mm_1 = (mm_hi_1 << 32) | mm_lo_1;
        uint64_t mm_2 = (mm_hi_2 << 32) | mm_lo_2;
        two_bits = (__uint128_t)mm_1 | ((__uint128_t)mm_2 << 64);
#endif
    }

    if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
        static const __m256i LUT_ACTG = make_lut_actg();
        __m256i mask_two_bits = _mm256_set1_epi8(0b110);
        __m256i mask_upper    = _mm256_set1_epi8((int8_t)0b11011111u);
        __m256i uv1 = _mm256_and_si256(v1, mask_upper);
        __m256i uv2 = _mm256_and_si256(v2, mask_upper);
        uint64_t mask_actg = movemask_64(
            _mm256_cmpeq_epi8(_mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v1, mask_two_bits)), uv1),
            _mm256_cmpeq_epi8(_mm256_shuffle_epi8(LUT_ACTG, _mm256_and_si256(v2, mask_two_bits)), uv2)
        );
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG))
            is_dna = mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG))
            mask_non_actg = ~mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_N)) {
            alignas(32) const char n[32] = {'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
                                            'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'};
            __m256i vn = _mm256_load_si256((const __m256i*)n);
            mask_n = u8_mask_avx2(v1, v2, vn);
        }
    }
}

} // anonymous namespace

template <uint64_t CONFIG>
inline FastaBitmask extract_fasta_bitmask(const uint8_t* buf) {
    using namespace advanced;
    __m256i v1 = _mm256_loadu_si256((const __m256i*)buf);
    __m256i v2 = _mm256_loadu_si256((const __m256i*)(buf + 32));

    alignas(32) const char gt[32] = {'>','>','>','>','>','>','>','>','>','>','>','>','>','>','>','>',
                                     '>','>','>','>','>','>','>','>','>','>','>','>','>','>','>','>'};
    alignas(32) const char lf[32] = {'\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n',
                                     '\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n'};
    __m256i vgt = _mm256_load_si256((const __m256i*)gt);
    __m256i vlf = _mm256_load_si256((const __m256i*)lf);

    FastaBitmask m{};
    m.is_dna = ~0ull;
    m.open_bracket = u8_mask_avx2(v1, v2, vgt);
    m.line_feeds   = u8_mask_avx2(v1, v2, vlf);
    bitpack_dna<CONFIG>(v1, v2, m.is_dna, m.two_bits, m.high_bit, m.low_bit, m.mask_non_actg, m.mask_n);
    return m;
}

template <uint64_t CONFIG>
inline FastqBitmask extract_fastq_bitmask(const uint8_t* buf) {
    using namespace advanced;
    __m256i v1 = _mm256_loadu_si256((const __m256i*)buf);
    __m256i v2 = _mm256_loadu_si256((const __m256i*)(buf + 32));

    alignas(32) const char lf[32] = {'\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n',
                                     '\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n','\n'};
    __m256i vlf = _mm256_load_si256((const __m256i*)lf);

    FastqBitmask m{};
    m.is_dna = ~0ull;
    m.line_feeds = u8_mask_avx2(v1, v2, vlf);
    bitpack_dna<CONFIG>(v1, v2, m.is_dna, m.two_bits, m.high_bit, m.low_bit, m.mask_non_actg, m.mask_n);
    return m;
}

} // namespace helicase::simd

// ─── SSSE3 path ──────────────────────────────────────────────────────────────
#elif defined(__SSSE3__)
#include <tmmintrin.h>

namespace helicase::simd {

namespace {

inline uint64_t movemask_64(__m128i v1, __m128i v2, __m128i v3, __m128i v4) {
    return ((uint64_t)(uint16_t)_mm_movemask_epi8(v1))         |
           (((uint64_t)(uint16_t)_mm_movemask_epi8(v2)) << 16) |
           (((uint64_t)(uint16_t)_mm_movemask_epi8(v3)) << 32) |
           (((uint64_t)(uint16_t)_mm_movemask_epi8(v4)) << 48);
}

inline uint64_t u8_mask_sse(__m128i v1, __m128i v2, __m128i v3, __m128i v4, __m128i vc) {
    return movemask_64(_mm_cmpeq_epi8(v1, vc), _mm_cmpeq_epi8(v2, vc),
                       _mm_cmpeq_epi8(v3, vc), _mm_cmpeq_epi8(v4, vc));
}

template <uint64_t CONFIG>
inline void bitpack_dna(__m128i v1, __m128i v2, __m128i v3, __m128i v4,
                        uint64_t& is_dna, __uint128_t& two_bits,
                        uint64_t& high_bit, uint64_t& low_bit,
                        uint64_t& mask_non_actg, uint64_t& mask_n)
{
    using namespace advanced;

    if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR)) {
        uint64_t hi1 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v1, 5));
        uint64_t lo1 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v1, 6));
        uint64_t hi2 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v2, 5));
        uint64_t lo2 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v2, 6));
        uint64_t hi3 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v3, 5));
        uint64_t lo3 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v3, 6));
        uint64_t hi4 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v4, 5));
        uint64_t lo4 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v4, 6));
        high_bit = hi1 | (hi2 << 16) | (hi3 << 32) | (hi4 << 48);
        low_bit  = lo1 | (lo2 << 16) | (lo3 << 32) | (lo4 << 48);
    }

    if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_PACKED)) {
#if defined(__BMI2__)
        uint64_t hi1 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v1, 5));
        uint64_t lo1 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v1, 6));
        uint64_t hi2 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v2, 5));
        uint64_t lo2 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v2, 6));
        uint64_t hi3 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v3, 5));
        uint64_t lo3 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v3, 6));
        uint64_t hi4 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v4, 5));
        uint64_t lo4 = (uint16_t)_mm_movemask_epi8(_mm_slli_epi16(v4, 6));
        uint64_t hi_12 = hi1 | (hi2 << 16), lo_12 = lo1 | (lo2 << 16);
        uint64_t hi_34 = hi3 | (hi4 << 16), lo_34 = lo3 | (lo4 << 16);
        uint64_t mm_1 = _pdep_u64(hi_12, 0xAAAAAAAAAAAAAAAAull) | _pdep_u64(lo_12, 0x5555555555555555ull);
        uint64_t mm_2 = _pdep_u64(hi_34, 0xAAAAAAAAAAAAAAAAull) | _pdep_u64(lo_34, 0x5555555555555555ull);
        two_bits = (__uint128_t)mm_1 | ((__uint128_t)mm_2 << 64);
#else
        // No-PDEP: unpack approach
        __m128i hi1 = _mm_slli_epi16(v1, 5), lo1 = _mm_slli_epi16(v1, 6);
        __m128i hi2 = _mm_slli_epi16(v2, 5), lo2 = _mm_slli_epi16(v2, 6);
        __m128i hi3 = _mm_slli_epi16(v3, 5), lo3 = _mm_slli_epi16(v3, 6);
        __m128i hi4 = _mm_slli_epi16(v4, 5), lo4 = _mm_slli_epi16(v4, 6);
        uint64_t mhi1 = (uint16_t)_mm_movemask_epi8(_mm_unpackhi_epi8(lo1, hi1));
        uint64_t mlo1 = (uint16_t)_mm_movemask_epi8(_mm_unpacklo_epi8(lo1, hi1));
        uint64_t mhi2 = (uint16_t)_mm_movemask_epi8(_mm_unpackhi_epi8(lo2, hi2));
        uint64_t mlo2 = (uint16_t)_mm_movemask_epi8(_mm_unpacklo_epi8(lo2, hi2));
        uint64_t mhi3 = (uint16_t)_mm_movemask_epi8(_mm_unpackhi_epi8(lo3, hi3));
        uint64_t mlo3 = (uint16_t)_mm_movemask_epi8(_mm_unpacklo_epi8(lo3, hi3));
        uint64_t mhi4 = (uint16_t)_mm_movemask_epi8(_mm_unpackhi_epi8(lo4, hi4));
        uint64_t mlo4 = (uint16_t)_mm_movemask_epi8(_mm_unpacklo_epi8(lo4, hi4));
        uint64_t mm_1 = mlo1 | (mhi1 << 16) | (mlo2 << 32) | (mhi2 << 48);
        uint64_t mm_2 = mlo3 | (mhi3 << 16) | (mlo4 << 32) | (mhi4 << 48);
        two_bits = (__uint128_t)mm_1 | ((__uint128_t)mm_2 << 64);
#endif
    }

    if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
        alignas(16) const char lut_s[16] = "A_C_T_G_________";
        __m128i LUT_ACTG = _mm_load_si128((const __m128i*)lut_s);
        __m128i m2b = _mm_set1_epi8(0b110);
        __m128i mup = _mm_set1_epi8((int8_t)0b11011111u);
        __m128i uv1 = _mm_and_si128(v1, mup), uv2 = _mm_and_si128(v2, mup);
        __m128i uv3 = _mm_and_si128(v3, mup), uv4 = _mm_and_si128(v4, mup);
        uint64_t mask_actg = movemask_64(
            _mm_cmpeq_epi8(_mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v1, m2b)), uv1),
            _mm_cmpeq_epi8(_mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v2, m2b)), uv2),
            _mm_cmpeq_epi8(_mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v3, m2b)), uv3),
            _mm_cmpeq_epi8(_mm_shuffle_epi8(LUT_ACTG, _mm_and_si128(v4, m2b)), uv4)
        );
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG))
            is_dna = mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG))
            mask_non_actg = ~mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_N)) {
            __m128i vn = _mm_set1_epi8('N');
            mask_n = u8_mask_sse(uv1, uv2, uv3, uv4, vn);
        }
    }
}

} // anonymous namespace

template <uint64_t CONFIG>
inline FastaBitmask extract_fasta_bitmask(const uint8_t* buf) {
    __m128i v1 = _mm_loadu_si128((const __m128i*)buf);
    __m128i v2 = _mm_loadu_si128((const __m128i*)(buf + 16));
    __m128i v3 = _mm_loadu_si128((const __m128i*)(buf + 32));
    __m128i v4 = _mm_loadu_si128((const __m128i*)(buf + 48));
    __m128i vgt = _mm_set1_epi8('>');
    __m128i vlf = _mm_set1_epi8('\n');
    FastaBitmask m{};
    m.is_dna = ~0ull;
    m.open_bracket = u8_mask_sse(v1, v2, v3, v4, vgt);
    m.line_feeds   = u8_mask_sse(v1, v2, v3, v4, vlf);
    bitpack_dna<CONFIG>(v1, v2, v3, v4, m.is_dna, m.two_bits, m.high_bit, m.low_bit, m.mask_non_actg, m.mask_n);
    return m;
}

template <uint64_t CONFIG>
inline FastqBitmask extract_fastq_bitmask(const uint8_t* buf) {
    __m128i v1 = _mm_loadu_si128((const __m128i*)buf);
    __m128i v2 = _mm_loadu_si128((const __m128i*)(buf + 16));
    __m128i v3 = _mm_loadu_si128((const __m128i*)(buf + 32));
    __m128i v4 = _mm_loadu_si128((const __m128i*)(buf + 48));
    __m128i vlf = _mm_set1_epi8('\n');
    FastqBitmask m{};
    m.is_dna = ~0ull;
    m.line_feeds = u8_mask_sse(v1, v2, v3, v4, vlf);
    bitpack_dna<CONFIG>(v1, v2, v3, v4, m.is_dna, m.two_bits, m.high_bit, m.low_bit, m.mask_non_actg, m.mask_n);
    return m;
}

} // namespace helicase::simd

// ─── Scalar fallback ─────────────────────────────────────────────────────────
#else

namespace helicase::simd {

static constexpr uint8_t UPPERCASE = 0b11011111u;
static constexpr uint8_t LUT_ACTG[16] = {
    'A','_','C','_','T','_','G','_','_','_','_','_','_','_','_','_'
};

template <uint64_t CONFIG>
inline FastaBitmask extract_fasta_bitmask(const uint8_t* buf) {
    using namespace advanced;
    FastaBitmask m{};
    m.is_dna = ~0ull;
    uint64_t mask_actg = 0;

    for (int i = 0; i < 64; ++i) {
        uint8_t x = buf[i];
        uint64_t bit = 1ull << i;
        m.open_bracket |= (x == '>') ? bit : 0;
        m.line_feeds   |= (x == '\n') ? bit : 0;

        if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR)) {
            m.high_bit |= (uint64_t)((x >> 2) & 1) << i;
            m.low_bit  |= (uint64_t)((x >> 1) & 1) << i;
        }
        if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_PACKED)) {
            m.two_bits |= (__uint128_t)((x >> 1) & 0b11) << (2 * i);
        }
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
            uint64_t is_actg = ((x & UPPERCASE) == LUT_ACTG[x & 0b110]) ? bit : 0;
            mask_actg |= is_actg;
        }
    }
    if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG))
            m.is_dna = mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG))
            m.mask_non_actg = ~mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_N)) {
            for (int i = 0; i < 64; ++i)
                m.mask_n |= (buf[i] == 'N') ? (1ull << i) : 0;
        }
    }
    return m;
}

template <uint64_t CONFIG>
inline FastqBitmask extract_fastq_bitmask(const uint8_t* buf) {
    using namespace advanced;
    FastqBitmask m{};
    m.is_dna = ~0ull;
    uint64_t mask_actg = 0;

    for (int i = 0; i < 64; ++i) {
        uint8_t x = buf[i];
        uint64_t bit = 1ull << i;
        m.line_feeds |= (x == '\n') ? bit : 0;

        if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR)) {
            m.high_bit |= (uint64_t)((x >> 2) & 1) << i;
            m.low_bit  |= (uint64_t)((x >> 1) & 1) << i;
        }
        if constexpr (flag_is_set(CONFIG, COMPUTE_DNA_PACKED)) {
            m.two_bits |= (__uint128_t)((x >> 1) & 0b11) << (2 * i);
        }
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
            uint64_t is_actg = ((x & UPPERCASE) == LUT_ACTG[x & 0b110]) ? bit : 0;
            mask_actg |= is_actg;
        }
    }
    if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG | COMPUTE_MASK_NON_ACTG)) {
        if constexpr (flag_is_set(CONFIG, SPLIT_NON_ACTG))
            m.is_dna = mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_NON_ACTG))
            m.mask_non_actg = ~mask_actg;
        if constexpr (flag_is_set(CONFIG, COMPUTE_MASK_N)) {
            for (int i = 0; i < 64; ++i)
                m.mask_n |= (buf[i] == 'N') ? (1ull << i) : 0;
        }
    }
    return m;
}

} // namespace helicase::simd

#endif // SIMD selection
