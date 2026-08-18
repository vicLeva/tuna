#pragma once

// On-disk / in-memory superkmer format:
//   [hdr_t len_bases][ceil(len/4) packed bytes]
//
// hdr_t is uint8_t when max superkmer length fits in 8 bits (2k-m ≤ 255),
// uint16_t otherwise.  The type is a compile-time constant deduced from k and m,
// so the header is as compact as possible without ever truncating a value.
//
//   Max superkmer length : 2k − m  (two k-mers sharing the boundary minimizer)
// For k ≤ 135 with m ≥ 15 (the common genomics range), hdr_t = uint8_t.
// For larger k, hdr_t = uint16_t.
//
// Packed encoding: base i is at bits 7-2*(i%4) of byte i/4, MSB-first.

#include <fstream>
#include <limits>
#include <mutex>
#include <string>
#include <cstdint>
#include <cstring>
#include <type_traits>

// Returns the path to the superkmer file for partition p under work_dir.
inline std::string partition_path(const std::string& work_dir, size_t p)
{
    return work_dir + "hash_" + std::to_string(p) + ".superkmers";
}

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef MAP_POPULATE
#  define MAP_POPULATE 0
#endif

// Deduce the header integer type for a given (k, m) pair.
// uint8_t when 2k − m ≤ 255, uint16_t otherwise.
template <uint16_t k, uint16_t m>
using sk_hdr_t = std::conditional_t<(2u * k - m <= 255u), uint8_t, uint16_t>;

// Compatibility sentinel: rolling-hash routing needs no minimizer coordinate.
template <uint16_t k, uint16_t m>
static constexpr sk_hdr_t<k, m> sk_no_min = std::numeric_limits<sk_hdr_t<k, m>>::max();


// ─── Writer ───────────────────────────────────────────────────────────────────
//
// SuperkmerWriter uses a raw POD buffer instead of std::string to avoid
// zero-initialisation on resize.  std::string::resize zeroes new bytes even when
// immediately overwritten; this buffer grows by doubling but skips zero-fill.

template <uint16_t k, uint16_t m>
struct SuperkmerWriter
{
    using hdr_t = sk_hdr_t<k, m>;              // superkmer header type (local alias)
    static constexpr size_t HDR_BYTES = sizeof(hdr_t);
    static constexpr size_t MAX_RECORD_BYTES =
        HDR_BYTES + (static_cast<size_t>(2u * k - m) + 3u) / 4u;

    char*  raw_  = nullptr;  // raw buffer pointer
    size_t sz_   = 0;        // used bytes
    size_t cap_  = 0;        // allocated bytes
    // flush_threshold is set per-writer based on n_parts so that total writer
    // memory across all partitions stays bounded to ~64 MB per thread:
    //   max(4 KB, 64 MB / n_parts)
    // Default (512 KB) is used when n_parts is small (≤128).
    size_t flush_threshold;

    explicit SuperkmerWriter(size_t flush_thresh = 512u << 10)
        // flush is checked after append, so retain one full-record of slack.
        // Otherwise every active writer eventually grows and doubles just
        // before the caller clears it.
        : cap_(flush_thresh + MAX_RECORD_BYTES), flush_threshold(flush_thresh)
    {
        raw_ = static_cast<char*>(::operator new(cap_));
    }

    // Copy constructor: used by vector(n, template) — allocates a fresh buffer
    // of the same capacity (preserving the pre-allocation invariant).
    SuperkmerWriter(const SuperkmerWriter& o)
        : sz_(o.sz_), cap_(o.cap_), flush_threshold(o.flush_threshold)
    {
        raw_ = static_cast<char*>(::operator new(cap_));
        if (sz_) std::memcpy(raw_, o.raw_, sz_);
    }

    SuperkmerWriter(SuperkmerWriter&& o) noexcept
        : raw_(o.raw_), sz_(o.sz_), cap_(o.cap_), flush_threshold(o.flush_threshold)
    { o.raw_ = nullptr; o.sz_ = 0; o.cap_ = 0; }

    SuperkmerWriter& operator=(const SuperkmerWriter&) = delete;
    SuperkmerWriter& operator=(SuperkmerWriter&&)      = delete;

    ~SuperkmerWriter() { ::operator delete(raw_); }

    const char* data()  const noexcept { return raw_; }
    size_t      size()  const noexcept { return sz_; }
    bool        empty() const noexcept { return sz_ == 0; }
    void        clear()       noexcept { sz_ = 0; }

    // Grow the buffer if sz_ + extra would exceed cap_.
    [[gnu::cold, gnu::noinline]]
    void grow(size_t extra)
    {
        const size_t need = sz_ + extra;
        cap_ = std::max(need, cap_ * 2);
        char* n = static_cast<char*>(::operator new(cap_));
        std::memcpy(n, raw_, sz_);
        ::operator delete(raw_);
        raw_ = n;
    }

    // Inline fast-path: capacity check + advance sz_.
    char* reserve_inline(size_t extra) noexcept
    {
        if (__builtin_expect(sz_ + extra > cap_, 0)) grow(extra);
        char* dst = raw_ + sz_;
        sz_ += extra;
        return dst;
    }

    // Serialise one superkmer from ASCII DNA (ACGT, any case).
    void append(const char* data, hdr_t len, hdr_t min_pos)
    {
        (void)min_pos;
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst, &len, sizeof(hdr_t));
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + HDR_BYTES);

        size_t i = 0;
        for (; i + 4 <= static_cast<size_t>(len); i += 4) {
            const uint8_t b0 = ((uint8_t(data[i  ]) >> 2) ^ (uint8_t(data[i  ]) >> 1)) & 3u;
            const uint8_t b1 = ((uint8_t(data[i+1]) >> 2) ^ (uint8_t(data[i+1]) >> 1)) & 3u;
            const uint8_t b2 = ((uint8_t(data[i+2]) >> 2) ^ (uint8_t(data[i+2]) >> 1)) & 3u;
            const uint8_t b3 = ((uint8_t(data[i+3]) >> 2) ^ (uint8_t(data[i+3]) >> 1)) & 3u;
            packed[i >> 2] = static_cast<uint8_t>((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
        }
        if (i < static_cast<size_t>(len)) {
            uint8_t tail = 0;
            for (size_t j = i; j < static_cast<size_t>(len); ++j) {
                const uint8_t b = ((uint8_t(data[j]) >> 2) ^ (uint8_t(data[j]) >> 1)) & 3u;
                tail |= static_cast<uint8_t>(b << (6u - 2u * (j - i)));
            }
            packed[i >> 2] = tail;
        }
    }

    // Serialise one superkmer from pre-encoded kache bytes (A=0,C=1,G=2,T=3).
    void append_kache(const uint8_t* kdata, hdr_t len, hdr_t min_pos)
    {
        (void)min_pos;
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst, &len, sizeof(hdr_t));
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + HDR_BYTES);

        size_t i = 0;
        for (; i + 4 <= static_cast<size_t>(len); i += 4)
            packed[i >> 2] = static_cast<uint8_t>(
                (kdata[i] << 6) | (kdata[i+1] << 4) | (kdata[i+2] << 2) | kdata[i+3]);
        if (i < static_cast<size_t>(len)) {
            uint8_t tail = 0;
            for (size_t j = i; j < static_cast<size_t>(len); ++j)
                tail |= static_cast<uint8_t>(kdata[j] << (6u - 2u * (j - i)));
            packed[i >> 2] = tail;
        }
    }

    // Serialise a base-aligned slice from an already packed sequence. `packed`
    // must have one zero sentinel byte after its final data byte.
    void append_packed(const uint8_t* packed, size_t base_offset, hdr_t len, hdr_t min_pos)
    {
        (void)min_pos;
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst, &len, sizeof(hdr_t));
        auto* out = reinterpret_cast<uint8_t*>(dst + HDR_BYTES);

        const auto* src = packed + (base_offset >> 2);
        const unsigned shift = static_cast<unsigned>((base_offset & 3u) * 2u);
        if (shift == 0) {
            std::memcpy(out, src, packed_bytes);
        } else {
            static constexpr uint64_t LOW_MASKS[4] = {
                0,
                0x3f3f3f3f3f3f3f3fULL,
                0x0f0f0f0f0f0f0f0fULL,
                0x0303030303030303ULL
            };
            const uint64_t low_mask = LOW_MASKS[shift >> 1];
            const uint64_t high_mask = ~low_mask;

            size_t i = 0;
            for (; i + sizeof(uint64_t) <= packed_bytes; i += sizeof(uint64_t)) {
                uint64_t lo, hi;
                std::memcpy(&lo, src + i,     sizeof(lo));
                std::memcpy(&hi, src + i + 1, sizeof(hi));
                const uint64_t shifted =
                    ((lo & low_mask) << shift) |
                    ((hi & high_mask) >> (8u - shift));
                std::memcpy(out + i, &shifted, sizeof(shifted));
            }
            for (; i < packed_bytes; ++i)
                out[i] = static_cast<uint8_t>(
                    (src[i] << shift) | (src[i + 1] >> (8u - shift)));
        }

        const unsigned tail_bases = static_cast<unsigned>(len & 3u);
        if (tail_bases != 0)
            out[packed_bytes - 1] &= static_cast<uint8_t>(0xffu << (8u - 2u * tail_bases));
    }

    bool needs_flush() const noexcept { return sz_ >= flush_threshold; }

    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(raw_, static_cast<std::streamsize>(sz_));
        sz_ = 0;
    }

    void flush_to_mem(std::string& dst, std::mutex& mtx)
    {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        dst.append(raw_, sz_);
        sz_ = 0;
    }
};


// ─── Reader ───────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
struct SuperkmerReader
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = sizeof(hdr_t);

    explicit SuperkmerReader(const std::string& path)
    {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return;

        struct stat sb;
        if (fstat(fd_, &sb) < 0) return;
        size_ = static_cast<size_t>(sb.st_size);
        if (size_ == 0) {
            ok_ = true;
            return;
        }

        void* p = mmap(nullptr, size_, PROT_READ,
                       MAP_PRIVATE | MAP_POPULATE, fd_, 0);
        if (p == MAP_FAILED) return;
        map_ = static_cast<const char*>(p);

        madvise(const_cast<char*>(map_), size_,
                MADV_SEQUENTIAL | MADV_WILLNEED);

        cur_ = map_;
        end_ = map_ + size_;
        ok_ = true;
    }

    ~SuperkmerReader()
    {
        if (map_) munmap(const_cast<char*>(map_), size_);
        if (fd_ >= 0) close(fd_);
    }

    SuperkmerReader(const SuperkmerReader&)            = delete;
    SuperkmerReader& operator=(const SuperkmerReader&) = delete;

    bool next()
    {
        if (cur_ == nullptr) return false;
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        record_ = cur_;
        hdr_t len;
        std::memcpy(&len, cur_, sizeof(hdr_t));
        if (len == 0) return false;
        cur_ += HDR_BYTES;

        const size_t packed_bytes = (static_cast<size_t>(len) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;

        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len;
        cur_ += packed_bytes;
        record_size_ = HDR_BYTES + packed_bytes;
        return true;
    }

    void           reset()       noexcept { cur_ = map_; }
    const uint8_t* packed_data() const { return ptr_; }
    const char*    record_data() const { return record_; }
    size_t         record_size() const { return record_size_; }
    size_t         size()        const { return len_; }
    hdr_t          min_pos()     const { return sk_no_min<k, m>; }
    size_t         file_size()   const { return size_; }
    size_t         data_size()   const { return size_; }
    bool           ok()          const { return ok_; }

private:
    int            fd_      = -1;      // file descriptor
    size_t         size_    = 0;       // file size in bytes
    const char*    map_     = nullptr; // mmap base pointer
    const char*    cur_     = nullptr; // read cursor
    const char*    end_     = nullptr; // one past last byte
    bool           ok_      = false;
    const uint8_t* ptr_     = nullptr; // packed data of current superkmer
    size_t         len_     = 0;       // length of current superkmer (bases)
    const char*    record_  = nullptr; // start of current encoded record
    size_t         record_size_ = 0;   // bytes in current encoded record
};


// ─── In-memory reader (streaming mode) ───────────────────────────────────────

template <uint16_t k, uint16_t m>
struct MemoryReader
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = sizeof(hdr_t);

    MemoryReader() = default;
    explicit MemoryReader(const std::string& data) noexcept
        : begin_(data.data()), cur_(data.data()), end_(data.data() + data.size()) {}

    bool next() noexcept
    {
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        record_ = cur_;
        hdr_t len;
        std::memcpy(&len, cur_, sizeof(hdr_t));
        if (len == 0) return false;
        cur_ += HDR_BYTES;
        const size_t packed_bytes = (static_cast<size_t>(len) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;
        ptr_  = reinterpret_cast<const uint8_t*>(cur_);
        len_  = len;
        cur_ += packed_bytes;
        record_size_ = HDR_BYTES + packed_bytes;
        return true;
    }

    void           reset()       noexcept { cur_ = begin_; }
    const uint8_t* packed_data() const noexcept { return ptr_; }
    const char*    record_data() const noexcept { return record_; }
    size_t         record_size() const noexcept { return record_size_; }
    size_t         size()        const noexcept { return len_; }
    hdr_t          min_pos()     const noexcept { return sk_no_min<k, m>; }
    bool           ok()          const noexcept { return cur_ != nullptr; }
    size_t         data_size()   const noexcept {
        return static_cast<size_t>(end_ - begin_);
    }

private:
    const char*    begin_   = nullptr;
    const char*    cur_     = nullptr;
    const char*    end_     = nullptr;
    const uint8_t* ptr_     = nullptr;
    size_t         len_     = 0;
    const char*    record_  = nullptr;
    size_t         record_size_ = 0;
};
