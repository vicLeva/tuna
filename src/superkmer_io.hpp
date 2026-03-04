#pragma once

// On-disk superkmer format: [uint32_t len][len bytes of sequence]
//
// SuperkmerWriter  — per-thread per-bucket buffered write.
//   Accumulates records in a local string buffer; caller flushes to
//   the shared ofstream under its mutex when needs_flush() is true.
//
// SuperkmerReader  — sequential read of one partition file.
//   Calling next() loads the next record; data()/size()/operator[] expose it.

#include <fstream>
#include <string>
#include <mutex>
#include <cstdint>
#include <cstring>


// ─── Writer ───────────────────────────────────────────────────────────────────

struct SuperkmerWriter
{
    std::string buf;
    static constexpr size_t FLUSH_THRESHOLD = 1u << 16; // 64 KB

    // Serialise one superkmer record into the local buffer.
    void append(const char* data, uint32_t len)
    {
        const size_t off = buf.size();
        buf.resize(off + sizeof(uint32_t) + len);
        std::memcpy(buf.data() + off,                    &len,  sizeof(uint32_t));
        std::memcpy(buf.data() + off + sizeof(uint32_t), data,  len);
    }

    bool needs_flush() const { return buf.size() >= FLUSH_THRESHOLD; }

    // Flush to the shared file under its mutex; no-op if buffer is empty.
    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(buf.data(), buf.size());
        buf.clear();
    }
};


// ─── Reader ───────────────────────────────────────────────────────────────────

struct SuperkmerReader
{
    explicit SuperkmerReader(const std::string& path)
        : in_(path, std::ios::binary) {}

    // Advance to the next superkmer.  Returns false at EOF or read error.
    bool next()
    {
        uint32_t len;
        if (!in_.read(reinterpret_cast<char*>(&len), sizeof(uint32_t)))
            return false;
        buf_.resize(len);
        in_.read(buf_.data(), len);
        return true;
    }

    const char* data()           const { return buf_.data(); }
    size_t      size()           const { return buf_.size(); }
    char        operator[](size_t i) const { return buf_[i]; }

private:
    std::ifstream in_;
    std::string   buf_;
};
