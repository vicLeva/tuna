#pragma once
#include <cstdint>
#include <cstring>
#include <ostream>
#include <vector>

namespace helicase {

// ─── ColumnarDNA ────────────────────────────────────────────────────────────
// Stores DNA as two separate bit vectors (high_bits and low_bits).
// Encoding: A=(0,0), C=(0,1), T=(1,0), G=(1,1) using (high_bit, low_bit).
class ColumnarDNA {
    std::vector<uint64_t> high_bits_;
    std::vector<uint64_t> low_bits_;
    uint64_t cur_hi_ = 0;
    uint64_t cur_lo_ = 0;
    size_t len_ = 0;

public:
    ColumnarDNA() = default;
    static ColumnarDNA with_capacity(size_t cap) {
        ColumnarDNA d;
        d.high_bits_.reserve(cap);
        d.low_bits_.reserve(cap);
        return d;
    }

    size_t len() const noexcept { return len_; }
    bool empty() const noexcept { return len_ == 0; }
    size_t capacity() const noexcept { return high_bits_.capacity(); }

    void clear() noexcept {
        high_bits_.clear();
        low_bits_.clear();
        cur_hi_ = cur_lo_ = 0;
        len_ = 0;
    }

    // Append `size` bases from packed bit vectors (high, low).
    // Bits 0..size-1 of high and low are used.
    void append(uint64_t high, uint64_t low, size_t size) noexcept {
        if (size == 0) return;
        size_t rem = len_ % 64;
        // mask the input to `size` bits
        uint64_t mask = (size >= 64) ? ~0ull : (~0ull >> (64 - size));
        len_ += size;
        uint64_t hi = high & mask;
        uint64_t lo = low & mask;
        if (rem + size >= 64) {
            cur_hi_ |= high << rem;  // rem can be 0 → shift by 0 is fine
            cur_lo_ |= low << rem;
            high_bits_.push_back(cur_hi_);
            low_bits_.push_back(cur_lo_);
            // Guard: hi >> (64 - rem) is UB when rem == 0; result would be 0.
            cur_hi_ = (rem == 0) ? 0ull : (hi >> (64 - rem));
            cur_lo_ = (rem == 0) ? 0ull : (lo >> (64 - rem));
        } else {
            cur_hi_ |= hi << rem;
            cur_lo_ |= lo << rem;
        }
    }

    // Returns (stored words, current partial word).
    std::pair<const std::vector<uint64_t>&, uint64_t> high_bits() const noexcept {
        return {high_bits_, cur_hi_};
    }
    std::pair<const std::vector<uint64_t>&, uint64_t> low_bits() const noexcept {
        return {low_bits_, cur_lo_};
    }

    // Get the (high, low) bits at index i.
    std::pair<bool, bool> get(size_t i) const noexcept {
        size_t word = i / 64, bit = i % 64;
        if (word < high_bits_.size()) {
            return {(high_bits_[word] >> bit) & 1, (low_bits_[word] >> bit) & 1};
        }
        return {(cur_hi_ >> bit) & 1, (cur_lo_ >> bit) & 1};
    }

    char get_char(size_t i) const noexcept {
        auto [h, l] = get(i);
        if (!h && !l) return 'A';
        if (!h && l)  return 'C';
        if (h  && !l) return 'T';
        return 'G';
    }

    friend std::ostream& operator<<(std::ostream& os, const ColumnarDNA& d) {
        for (size_t i = 0; i < d.len(); ++i) os << d.get_char(i);
        return os;
    }
};

// ─── PackedDNA ──────────────────────────────────────────────────────────────
// 2-bit packed DNA: A=0, C=1, T=2, G=3 (bits 0-1 per base).
// Stored in little-endian order (base 0 at bits 0-1 of first word).
class PackedDNA {
    using T = __uint128_t;
    static constexpr size_t BITS_PER_BLOCK = 128;
    std::vector<T> bits_;
    T cur_ = 0;
    size_t num_bits_ = 0;  // in bits (not bases)

public:
    PackedDNA() = default;
    static PackedDNA with_capacity(size_t cap) {
        PackedDNA d;
        d.bits_.reserve(cap);
        return d;
    }

    size_t len() const noexcept { return num_bits_ / 2; }
    bool empty() const noexcept { return num_bits_ == 0; }
    size_t capacity() const noexcept { return bits_.capacity(); }

    void clear() noexcept {
        bits_.clear();
        cur_ = 0;
        num_bits_ = 0;
    }

    // Append `num_bits` bits from `packed` (2 bits per base).
    void append(T packed, size_t num_bits) noexcept {
        if (num_bits == 0) return;
        size_t rem = num_bits_ % BITS_PER_BLOCK;
        T mask = (num_bits >= BITS_PER_BLOCK) ? ~(T)0 : (~(T)0 >> (BITS_PER_BLOCK - num_bits));
        num_bits_ += num_bits;
        T x = packed & mask;
        if (rem + num_bits >= BITS_PER_BLOCK) {
            cur_ |= packed << rem;  // rem==0: shift by 0 is fine
            bits_.push_back(cur_);
            cur_ = (rem == 0) ? (T)0 : (x >> (BITS_PER_BLOCK - rem));
        } else {
            cur_ |= x << rem;
        }
    }

    uint8_t get(size_t i) const noexcept {
        size_t word = i / 64, bit = i % 64;  // 64 bases per 128-bit word
        if (word < bits_.size()) {
            return (uint8_t)((bits_[word] >> (2 * bit)) & 0b11);
        }
        return (uint8_t)((cur_ >> (2 * bit)) & 0b11);
    }

    char get_char(size_t i) const noexcept {
        static constexpr char LUT[4] = {'A', 'C', 'T', 'G'};
        return LUT[get(i)];
    }

    // Raw access for downstream consumers.
    std::pair<const std::vector<T>&, T> bits() const noexcept {
        return {bits_, cur_};
    }

    friend std::ostream& operator<<(std::ostream& os, const PackedDNA& d) {
        for (size_t i = 0; i < d.len(); ++i) os << d.get_char(i);
        return os;
    }
};

// ─── BitMask ─────────────────────────────────────────────────────────────────
// Single-bit-per-position mask stored as packed uint64_t words.
class BitMask {
    std::vector<uint64_t> bits_;
    uint64_t cur_ = 0;
    size_t len_ = 0;

public:
    BitMask() = default;
    static BitMask with_capacity(size_t cap) {
        BitMask m;
        m.bits_.reserve(cap);
        return m;
    }

    size_t len() const noexcept { return len_; }
    bool empty() const noexcept { return len_ == 0; }
    size_t capacity() const noexcept { return bits_.capacity(); }

    void clear() noexcept {
        bits_.clear();
        cur_ = 0;
        len_ = 0;
    }

    void append(uint64_t x, size_t size) noexcept {
        if (size == 0) return;
        size_t rem = len_ % 64;
        uint64_t mask = (size >= 64) ? ~0ull : (~0ull >> (64 - size));
        len_ += size;
        uint64_t y = x & mask;
        if (rem + size >= 64) {
            cur_ |= x << rem;
            bits_.push_back(cur_);
            cur_ = (rem == 0) ? 0ull : (y >> (64 - rem));
        } else {
            cur_ |= y << rem;
        }
    }

    bool get(size_t i) const noexcept {
        size_t word = i / 64, bit = i % 64;
        if (word < bits_.size()) return (bits_[word] >> bit) & 1;
        return (cur_ >> bit) & 1;
    }

    std::pair<const std::vector<uint64_t>&, uint64_t> bits() const noexcept {
        return {bits_, cur_};
    }
};

} // namespace helicase
