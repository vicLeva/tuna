#pragma once

// On-disk superkmer format: [uint16_t len_bases][ceil(len/4) packed bytes]
//
// Each superkmer is stored as a 2-byte length prefix (number of bases, native
// endian) followed by the bases packed 4 per byte in kache-hash encoding
// (A=0, C=1, G=2, T=3), MSB-first (base i is at bits 7-2*(i%4) of byte i/4).
//
// This is 4x smaller than ASCII and avoids the ASCII↔base conversion in
// Phase 2: the unpacked value is a DNA::Base (kache-hash encoding) directly,
// consumable by Kmer_Window::advance(DNA::Base) without remapping.
//
// SuperkmerWriter — per-thread per-bucket buffered write.
//   Converts ASCII → packed on the fly; flushes to the shared ofstream under
//   its mutex when needs_flush() is true.
//
// SuperkmerReader — zero-copy sequential reader backed by mmap (Linux/POSIX).
//   Construction maps the entire file into the virtual address space.
//   next() reads the length prefix and advances the cursor; packed_data() and
//   size() return a pointer into the mapped region and the base count.

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
    static constexpr size_t FLUSH_THRESHOLD = 512u << 10; // 512 KB

    // Serialise one superkmer.
    // `data` is ASCII DNA (ACGT, any case); `len` is the number of bases.
    void append(const char* data, uint32_t len)
    {
        const size_t packed_bytes = (len + 3u) / 4u;
        const size_t off = buf.size();
        buf.resize(off + sizeof(uint16_t) + packed_bytes);

        // Length prefix.
        const uint16_t len16 = static_cast<uint16_t>(len);
        std::memcpy(buf.data() + off, &len16, sizeof(uint16_t));

        // Pack 4 bases per byte, kache encoding: ((c>>2)^(c>>1))&3 = A=0,C=1,G=2,T=3.
        uint8_t* packed = reinterpret_cast<uint8_t*>(buf.data() + off + sizeof(uint16_t));
        std::memset(packed, 0, packed_bytes);
        for (uint32_t i = 0; i < len; ++i) {
            const uint8_t b = ((uint8_t(data[i]) >> 2) ^ (uint8_t(data[i]) >> 1)) & 3u;
            packed[i >> 2] |= static_cast<uint8_t>(b << (6u - 2u * (i & 3u)));
        }
    }

    bool needs_flush() const { return buf.size() >= FLUSH_THRESHOLD; }

    // Flush to the shared file under its mutex; no-op if buffer is empty.
    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(buf.data(), static_cast<std::streamsize>(buf.size()));
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
        if (cur_ + static_cast<ptrdiff_t>(sizeof(uint16_t)) > end_) return false;
        uint16_t len16;
        std::memcpy(&len16, cur_, sizeof(uint16_t));
        if (len16 == 0) return false;
        cur_ += sizeof(uint16_t);

        const size_t packed_bytes = (static_cast<size_t>(len16) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;

        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len16;
        cur_ += packed_bytes;
        return true;
    }

    // Pointer to the packed bases of the current superkmer (kache encoding,
    // 4 bases/byte MSB-first).  Valid until the next call to next().
    const uint8_t* packed_data() const { return ptr_; }

    // Number of bases in the current superkmer.
    size_t size() const { return len_; }

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
    size_t      len_  = 0;
};
