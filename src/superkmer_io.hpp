#pragma once

// On-disk superkmer format: [sequence bytes]['\0']
//
// Each superkmer is stored as its raw DNA characters followed by a null byte.
// No length prefix — the null sentinel is located with memchr, which the CPU
// can execute in a single SIMD pass (glibc's AVX2 memchr).
//
// SuperkmerWriter — per-thread per-bucket buffered write.
//   Appends [data\0] to a local buffer; caller flushes to the shared ofstream
//   under its mutex when needs_flush() is true.
//
// SuperkmerReader — zero-copy sequential reader backed by mmap (Linux/POSIX).
//   Construction maps the entire file into the virtual address space.
//   next() advances the cursor with a single memchr; data()/size()/operator[]
//   return pointers directly into the mapped region — no copies, no allocs.

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
#  define MAP_POPULATE 0  // not available outside Linux; no-op
#endif


// ─── Writer ───────────────────────────────────────────────────────────────────

struct SuperkmerWriter
{
    std::string buf;
    static constexpr size_t FLUSH_THRESHOLD = 512u << 10; // 512 KB

    // Serialise one superkmer record as [data\0].
    void append(const char* data, uint32_t len)
    {
        const size_t off = buf.size();
        buf.resize(off + len + 1);
        std::memcpy(buf.data() + off, data, len);
        buf[off + len] = '\0';
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

        // Hint: we'll scan the file exactly once, front to back.
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

    // Non-copyable (owns the mapping).
    SuperkmerReader(const SuperkmerReader&)            = delete;
    SuperkmerReader& operator=(const SuperkmerReader&) = delete;

    // Advance to the next superkmer.  Returns false at EOF.
    bool next()
    {
        if (cur_ >= end_) return false;
        ptr_ = cur_;
        const char* nul = static_cast<const char*>(
            std::memchr(cur_, '\0', static_cast<size_t>(end_ - cur_)));
        if (!nul) return false;
        len_ = static_cast<size_t>(nul - cur_);
        cur_ = nul + 1;
        return len_ > 0;  // skip any spurious double-nulls
    }

    const char* data()           const { return ptr_; }
    size_t      size()           const { return len_; }
    char        operator[](size_t i) const { return ptr_[i]; }

    // Total mapped bytes — used by count.hpp to pre-size the hash table.
    size_t file_size()           const { return size_; }

    bool ok()                    const { return map_ != nullptr; }

private:
    int         fd_  = -1;
    size_t      size_ = 0;
    const char* map_  = nullptr;
    const char* cur_  = nullptr;
    const char* end_  = nullptr;
    const char* ptr_  = nullptr;
    size_t      len_  = 0;
};
