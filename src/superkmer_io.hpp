#pragma once

// On-disk superkmer format: [uint8_t len_bases][uint8_t min_pos][ceil(len/4) packed bytes]
//
// Each superkmer is stored as two header bytes followed by the bases packed 4
// per byte in kache-hash encoding (A=0, C=1, G=2, T=3), MSB-first.
//
//   len_bases — number of bases (max 255; superkmers are at most 2k−m bases).
//   min_pos   — 0-indexed start of the minimizer m-mer within the superkmer.
//               Stored so Phase 2 can compute ntHash(minimizer) in O(m) instead
//               of O(k) via MinimizerWindow::reset().
//
// Packed encoding: base i is at bits 7-2*(i%4) of byte i/4.

#include <fstream>
#include <mutex>
#include <string>
#include <cstdint>
#include <cstring>

// Returns the path to the superkmer file for partition p under work_dir.
// Used in pipeline setup, counting, and cleanup — centralises the naming convention.
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


// ─── Writer ───────────────────────────────────────────────────────────────────
//
// SuperkmerWriter uses a raw POD buffer instead of std::string to avoid
// zero-initialisation on resize.  std::string::resize zeroes new bytes even when
// immediately overwritten; this buffer grows by doubling but skips zero-fill.

struct SuperkmerWriter
{
    char*  raw_  = nullptr;
    size_t sz_   = 0;
    size_t cap_  = 0;
    // flush_threshold is set per-writer based on n_parts so that total writer
    // memory across all partitions stays bounded to ~64 MB per thread:
    //   max(4 KB, 64 MB / n_parts)
    // Default (512 KB) is used when n_parts is small (≤128).
    size_t flush_threshold;

    explicit SuperkmerWriter(size_t flush_thresh = 512u << 10)
        : cap_(flush_thresh), flush_threshold(flush_thresh)
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
    // Only called on the slow path (amortised O(1): doubles each time).
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
    // Returns pointer to the newly reserved region (caller fills it).
    char* reserve_inline(size_t extra) noexcept
    {
        if (__builtin_expect(sz_ + extra > cap_, 0)) grow(extra);
        char* dst = raw_ + sz_;
        sz_ += extra;
        return dst;
    }

    // Serialise one superkmer from ASCII DNA (ACGT, any case).
    // Kept for call sites that don't pre-encode; prefer append_kache when possible.
    void append(const char* data, uint8_t len, uint8_t min_pos)
    {
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(2u + packed_bytes);
        dst[0] = static_cast<char>(len);
        dst[1] = static_cast<char>(min_pos);
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + 2);

        // Direct 4-bases-per-byte writes: no memset, no read-modify-write in the
        // main loop → compiler auto-vectorizes (no loop-carried dependency).
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
    // Caller pre-converts the entire sequence ASCII→kache once; this path
    // skips all re-encoding — only shift+OR to pack 4 bases into 1 byte.
    void append_kache(const uint8_t* kdata, uint8_t len, uint8_t min_pos)
    {
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(2u + packed_bytes);
        dst[0] = static_cast<char>(len);
        dst[1] = static_cast<char>(min_pos);
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + 2);

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

    bool needs_flush() const noexcept { return sz_ >= flush_threshold; }

    // Flush to the shared file under its mutex; no-op if buffer is empty.
    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(raw_, static_cast<std::streamsize>(sz_));
        sz_ = 0;
    }

    // Flush to an in-memory string sink (streaming mode — avoids disk I/O).
    void flush_to_mem(std::string& dst, std::mutex& mtx)
    {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        dst.append(raw_, sz_);
        sz_ = 0;
    }
};


// ─── Reader ───────────────────────────────────────────────────────────────────

struct SuperkmerReader
{
    explicit SuperkmerReader(const std::string& path)
    {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return;

        struct stat sb;
        if (fstat(fd_, &sb) < 0 || sb.st_size == 0) return;
        size_ = static_cast<size_t>(sb.st_size);

        void* m = mmap(nullptr, size_, PROT_READ,
                       MAP_PRIVATE | MAP_POPULATE, fd_, 0);
        if (m == MAP_FAILED) return;
        map_ = static_cast<const char*>(m);

        madvise(const_cast<char*>(map_), size_,
                MADV_SEQUENTIAL | MADV_WILLNEED);

        cur_ = map_;
        end_ = map_ + size_;
    }

    ~SuperkmerReader()
    {
        if (map_) munmap(const_cast<char*>(map_), size_);
        if (fd_ >= 0) close(fd_);
    }

    SuperkmerReader(const SuperkmerReader&)            = delete;
    SuperkmerReader& operator=(const SuperkmerReader&) = delete;

    // Advance to the next superkmer.  Returns false at EOF.
    bool next()
    {
        if (cur_ + 2 > end_) return false;
        const uint8_t len8 = static_cast<uint8_t>(*cur_);
        if (len8 == 0) return false;
        min_pos_ = static_cast<uint8_t>(cur_[1]);
        cur_ += 2;

        const size_t packed_bytes = (static_cast<size_t>(len8) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;

        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len8;
        cur_ += packed_bytes;
        return true;
    }

    // Pointer to the packed bases of the current superkmer (kache encoding,
    // 4 bases/byte MSB-first).  Valid until the next call to next().
    const uint8_t* packed_data() const { return ptr_; }

    // Number of bases in the current superkmer.
    size_t size() const { return len_; }

    // 0-indexed start position of the minimizer m-mer within the current superkmer.
    uint8_t min_pos() const { return min_pos_; }

    // Total mapped bytes — available for diagnostics / capacity estimation.
    size_t file_size() const { return size_; }

    bool ok() const { return map_ != nullptr; }

private:
    int         fd_   = -1;
    size_t      size_ = 0;
    const char* map_  = nullptr;
    const char* cur_  = nullptr;
    const char* end_  = nullptr;
    const uint8_t* ptr_ = nullptr;
    size_t      len_     = 0;
    uint8_t     min_pos_ = 0;
};


// ─── In-memory reader (streaming mode) ───────────────────────────────────────
//
// Same interface as SuperkmerReader but backed by an existing std::string
// rather than an mmap'd file.  Used by the streaming pipeline to avoid the
// disk write + mmap round-trip between Phase 1 and Phase 2.

struct MemoryReader
{
    MemoryReader() = default;
    explicit MemoryReader(const std::string& data) noexcept
        : cur_(data.data()), end_(data.data() + data.size()) {}

    bool next() noexcept
    {
        if (cur_ + 2 > end_) return false;
        const uint8_t len8 = static_cast<uint8_t>(*cur_);
        if (len8 == 0) return false;
        min_pos_ = static_cast<uint8_t>(cur_[1]);
        cur_ += 2;
        const size_t packed_bytes = (static_cast<size_t>(len8) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;
        ptr_  = reinterpret_cast<const uint8_t*>(cur_);
        len_  = len8;
        cur_ += packed_bytes;
        return true;
    }

    const uint8_t* packed_data() const noexcept { return ptr_; }
    size_t         size()        const noexcept { return len_; }
    uint8_t        min_pos()     const noexcept { return min_pos_; }
    bool           ok()          const noexcept { return cur_ != nullptr; }

private:
    const char*    cur_     = nullptr;
    const char*    end_     = nullptr;
    const uint8_t* ptr_     = nullptr;
    size_t         len_     = 0;
    uint8_t        min_pos_ = 0;
};
