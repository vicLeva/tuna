#pragma once

// On-disk superkmer format: [uint8_t len_bases][uint8_t min_pos][ceil(len/4) packed bytes]
//
// Each superkmer is stored as two header bytes followed by the bases packed 4
// per byte in kache-hash encoding (A=0, C=1, G=2, T=3), MSB-first.
//
//   len_bases — number of bases (max 255; superkmers are at most 2k−l bases).
//   min_pos   — 0-indexed start position of the minimizer l-mer within the
//               superkmer.  Stored so Phase 2 can compute ntHash(minimizer)
//               in O(l) without running MinimizerWindow::reset() in O(k).
//
// Packed encoding: base i is at bits 7-2*(i%4) of byte i/4.
// This is 4x smaller than ASCII.  Phase 2 unpacks directly to DNA::Base
// (kache encoding) without any ASCII round-trip.
//
// SuperkmerWriter — per-thread per-bucket buffered write.
//   Converts ASCII → packed on the fly; flushes to the shared ofstream under
//   its mutex when needs_flush() is true.
//
// SuperkmerReader — zero-copy sequential reader backed by mmap (Linux/POSIX).
//   Construction maps the entire file into the virtual address space.
//   next() reads the two header bytes and advances the cursor; packed_data(),
//   size(), and min_pos() expose the current superkmer.

#include <fstream>
#include <mutex>
#include <string>
#include <cstdint>
#include <cstring>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef MAP_POPULATE
#  define MAP_POPULATE 0
#endif


// ─── Writer ───────────────────────────────────────────────────────────────────

struct SuperkmerWriter
{
    std::string buf;
    // flush_threshold is set per-writer based on n_parts so that total writer
    // memory across all partitions stays bounded to ~64 MB per thread:
    //   max(4 KB, 64 MB / n_parts)
    // Default (512 KB) is used when n_parts is small (≤128).
    const size_t flush_threshold;

    explicit SuperkmerWriter(size_t flush_thresh = 512u << 10)
        : flush_threshold(flush_thresh) {}

    // Serialise one superkmer.
    // `data` is ASCII DNA (ACGT, any case); `len` is the number of bases;
    // `min_pos` is the 0-indexed start of the minimizer l-mer within the superkmer.
    void append(const char* data, uint8_t len, uint8_t min_pos)
    {
        const size_t packed_bytes = (len + 3u) / 4u;
        const size_t off = buf.size();
        buf.resize(off + 2u + packed_bytes);

        // Two header bytes: length then minimizer position.
        buf[off]     = static_cast<char>(len);
        buf[off + 1] = static_cast<char>(min_pos);

        // Pack 4 bases per byte, kache encoding: ((c>>2)^(c>>1))&3 = A=0,C=1,G=2,T=3.
        uint8_t* packed = reinterpret_cast<uint8_t*>(buf.data() + off + 2u);
        std::memset(packed, 0, packed_bytes);
        for (uint8_t i = 0; i < len; ++i) {
            const uint8_t b = ((uint8_t(data[i]) >> 2) ^ (uint8_t(data[i]) >> 1)) & 3u;
            packed[i >> 2] |= static_cast<uint8_t>(b << (6u - 2u * (i & 3u)));
        }
    }

    bool needs_flush() const { return buf.size() >= flush_threshold; }

    // Flush to the shared file under its mutex; no-op if buffer is empty.
    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(buf.data(), static_cast<std::streamsize>(buf.size()));
        buf.clear();
    }

    // Flush to an in-memory string sink (streaming mode — avoids disk I/O).
    void flush_to_mem(std::string& dst, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        dst.append(buf);
        buf.clear();
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

    // 0-indexed start position of the minimizer l-mer within the current superkmer.
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
