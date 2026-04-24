#pragma once

// Protein superkmer on-disk / in-memory format:
//   [hdr_t len_aa][hdr_t min_pos][ceil(5*len/8) packed bytes]
//
// hdr_t: uint8_t when 2k-m ≤ 255, uint16_t otherwise (same policy as DNA).
//
// Packing: amino acid i occupies bits [5*i .. 5*i+4] of the bitstream,
// MSB-first within each byte (bit 7 = earliest).
// write5(buf, 5*i, val) / read5(buf, 5*i) implement the 5-bit field I/O.

#include "aa_encode.hpp"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <string>
#include <type_traits>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef MAP_POPULATE
#  define MAP_POPULATE 0
#endif

namespace prot {

// ── Header type ───────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
using psk_hdr_t = std::conditional_t<(2u * k - m <= 255u), uint8_t, uint16_t>;

template <uint16_t k, uint16_t m>
static constexpr psk_hdr_t<k, m> psk_no_min = std::numeric_limits<psk_hdr_t<k, m>>::max();


// ── 5-bit bitstream primitives ────────────────────────────────────────────────
//
// Both functions address a bit stream where bit 0 is the MSB of byte 0.
// A value at bit_pos occupies bits [bit_pos .. bit_pos+4].

inline void write5(uint8_t* buf, uint32_t bit_pos, uint8_t val5) noexcept {
    const uint32_t b      = bit_pos >> 3;
    const uint32_t off    = bit_pos & 7;
    const uint32_t hi     = (off <= 3) ? 5 : (8 - off);  // bits in byte b
    const uint32_t lo     = 5 - hi;                        // bits in byte b+1

    // Clear and write the hi most-significant bits of val5 into byte b.
    const uint32_t shift0 = 8 - off - hi;
    const uint8_t  mask0  = static_cast<uint8_t>(((1u << hi) - 1u) << shift0);
    buf[b] = (buf[b] & ~mask0) | static_cast<uint8_t>((val5 >> lo) << shift0);

    if (lo > 0) {
        // Write the lo least-significant bits of val5 into the top of byte b+1.
        const uint8_t mask1 = static_cast<uint8_t>(((1u << lo) - 1u) << (8 - lo));
        buf[b + 1] = (buf[b + 1] & ~mask1) | static_cast<uint8_t>(val5 << (8 - lo));
    }
}

inline uint8_t read5(const uint8_t* buf, uint32_t bit_pos) noexcept {
    const uint32_t b   = bit_pos >> 3;
    const uint32_t off = bit_pos & 7;
    const uint32_t hi  = (off <= 3) ? 5 : (8 - off);
    const uint32_t lo  = 5 - hi;
    const uint32_t shift0 = 8 - off - hi;
    uint8_t val = static_cast<uint8_t>((buf[b] >> shift0) & ((1u << hi) - 1u)) << lo;
    if (lo > 0)
        val |= static_cast<uint8_t>(buf[b + 1] >> (8 - lo)) & static_cast<uint8_t>((1u << lo) - 1u);
    return val;
}


// ── Writer ────────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct ProtSuperKmerWriter {
    using hdr_t = psk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    char*  raw_  = nullptr;
    size_t sz_   = 0;
    size_t cap_  = 0;
    size_t flush_threshold;

    explicit ProtSuperKmerWriter(size_t flush_thresh = 512u << 10)
        : cap_(flush_thresh), flush_threshold(flush_thresh)
    { raw_ = static_cast<char*>(::operator new(cap_)); }

    ProtSuperKmerWriter(const ProtSuperKmerWriter& o)
        : sz_(o.sz_), cap_(o.cap_), flush_threshold(o.flush_threshold)
    {
        raw_ = static_cast<char*>(::operator new(cap_));
        if (sz_) std::memcpy(raw_, o.raw_, sz_);
    }

    ProtSuperKmerWriter(ProtSuperKmerWriter&& o) noexcept
        : raw_(o.raw_), sz_(o.sz_), cap_(o.cap_), flush_threshold(o.flush_threshold)
    { o.raw_ = nullptr; o.sz_ = 0; o.cap_ = 0; }

    ProtSuperKmerWriter& operator=(const ProtSuperKmerWriter&) = delete;
    ProtSuperKmerWriter& operator=(ProtSuperKmerWriter&&)      = delete;

    ~ProtSuperKmerWriter() { ::operator delete(raw_); }

    const char* data()  const noexcept { return raw_; }
    size_t      size()  const noexcept { return sz_; }
    bool        empty() const noexcept { return sz_ == 0; }
    void        clear()       noexcept { sz_ = 0; }

    [[gnu::cold, gnu::noinline]]
    void grow(size_t extra) {
        const size_t need = sz_ + extra;
        cap_ = std::max(need, cap_ * 2);
        char* n = static_cast<char*>(::operator new(cap_));
        std::memcpy(n, raw_, sz_);
        ::operator delete(raw_);
        raw_ = n;
    }

    char* reserve_inline(size_t extra) noexcept {
        if (__builtin_expect(sz_ + extra > cap_, 0)) grow(extra);
        char* dst = raw_ + sz_;
        sz_ += extra;
        return dst;
    }

    // Serialise one superkmer from pre-encoded amino acids (5-bit values, 0-19).
    void append(const uint8_t* enc, hdr_t len, hdr_t min_pos) {
        const size_t packed_bytes = (static_cast<size_t>(len) * 5 + 7) / 8;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst,                 &len,     sizeof(hdr_t));
        std::memcpy(dst + sizeof(hdr_t), &min_pos, sizeof(hdr_t));
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + HDR_BYTES);
        // Zero the packed region first (write5 ORs into it).
        std::memset(packed, 0, packed_bytes);
        for (uint32_t i = 0; i < static_cast<uint32_t>(len); ++i)
            write5(packed, i * 5, enc[i]);
    }

    bool needs_flush() const noexcept { return sz_ >= flush_threshold; }

    void flush_to(std::ofstream& file, std::mutex& mtx) {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(raw_, static_cast<std::streamsize>(sz_));
        sz_ = 0;
    }

    void flush_to_mem(std::string& dst, std::mutex& mtx) {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        dst.append(raw_, sz_);
        sz_ = 0;
    }
};


// ── Disk reader (mmap) ────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct ProtSuperKmerReader {
    using hdr_t = psk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    explicit ProtSuperKmerReader(const std::string& path) {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return;
        struct stat sb;
        if (fstat(fd_, &sb) < 0 || sb.st_size == 0) return;
        size_ = static_cast<size_t>(sb.st_size);
        void* p = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd_, 0);
        if (p == MAP_FAILED) return;
        map_ = static_cast<const char*>(p);
        madvise(const_cast<char*>(map_), size_, MADV_SEQUENTIAL | MADV_WILLNEED);
        cur_ = map_;
        end_ = map_ + size_;
    }

    ~ProtSuperKmerReader() {
        if (map_) munmap(const_cast<char*>(map_), size_);
        if (fd_ >= 0) close(fd_);
    }

    ProtSuperKmerReader(const ProtSuperKmerReader&)            = delete;
    ProtSuperKmerReader& operator=(const ProtSuperKmerReader&) = delete;

    bool next() noexcept {
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        hdr_t len, mp;
        std::memcpy(&len, cur_,                sizeof(hdr_t));
        std::memcpy(&mp,  cur_ + sizeof(hdr_t), sizeof(hdr_t));
        if (len == 0) return false;
        min_pos_ = mp;
        cur_ += HDR_BYTES;
        const size_t pb = (static_cast<size_t>(len) * 5 + 7) / 8;
        if (cur_ + static_cast<ptrdiff_t>(pb) > end_) return false;
        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len;
        cur_ += pb;
        return true;
    }

    const uint8_t* packed_data() const noexcept { return ptr_; }
    size_t         size()        const noexcept { return len_; }
    hdr_t          min_pos()     const noexcept { return min_pos_; }
    bool           ok()          const noexcept { return map_ != nullptr; }

private:
    int            fd_      = -1;
    size_t         size_    = 0;
    const char*    map_     = nullptr;
    const char*    cur_     = nullptr;
    const char*    end_     = nullptr;
    const uint8_t* ptr_     = nullptr;
    size_t         len_     = 0;
    hdr_t          min_pos_ = 0;
};


// ── In-memory reader ──────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct ProtMemoryReader {
    using hdr_t = psk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    ProtMemoryReader() = default;
    explicit ProtMemoryReader(const std::string& data) noexcept
        : cur_(data.data()), end_(data.data() + data.size()) {}

    bool next() noexcept {
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        hdr_t len, mp;
        std::memcpy(&len, cur_,                sizeof(hdr_t));
        std::memcpy(&mp,  cur_ + sizeof(hdr_t), sizeof(hdr_t));
        if (len == 0) return false;
        min_pos_ = mp;
        cur_ += HDR_BYTES;
        const size_t pb = (static_cast<size_t>(len) * 5 + 7) / 8;
        if (cur_ + static_cast<ptrdiff_t>(pb) > end_) return false;
        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len;
        cur_ += pb;
        return true;
    }

    const uint8_t* packed_data() const noexcept { return ptr_; }
    size_t         size()        const noexcept { return len_; }
    hdr_t          min_pos()     const noexcept { return min_pos_; }
    bool           ok()          const noexcept { return cur_ != nullptr; }

private:
    const char*    cur_     = nullptr;
    const char*    end_     = nullptr;
    const uint8_t* ptr_     = nullptr;
    size_t         len_     = 0;
    hdr_t          min_pos_ = 0;
};


// ── Partition file path ───────────────────────────────────────────────────────

inline std::string prot_partition_path(const std::string& work_dir, size_t p)
{
    return work_dir + "prot_" + std::to_string(p) + ".superkmers";
}

} // namespace prot
