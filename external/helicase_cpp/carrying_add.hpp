#pragma once
#include <cstdint>
#ifdef __x86_64__
#  include <immintrin.h>
#endif

namespace helicase {

// Propagates carry across 64-byte block boundaries for header detection.
// Uses _addcarry_u64 on x86_64 for efficiency.
struct Carry {
#ifdef __x86_64__
    uint8_t c = 0;
    explicit constexpr Carry(bool v = false) noexcept : c(v ? 1 : 0) {}

    inline uint64_t add(uint64_t lhs, uint64_t rhs) noexcept {
        unsigned long long res;
        c = _addcarry_u64(c, lhs, rhs, &res);
        return (uint64_t)res;
    }
#else
    bool c = false;
    explicit constexpr Carry(bool v = false) noexcept : c(v) {}

    inline uint64_t add(uint64_t lhs, uint64_t rhs) noexcept {
        uint64_t a = lhs + rhs;
        bool c1 = a < lhs;
        uint64_t b = a + (uint64_t)c;
        bool c2 = b < a;
        c = c1 | c2;
        return b;
    }
#endif
};

} // namespace helicase
