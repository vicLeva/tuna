#pragma once
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef _WIN32
#  include <sys/mman.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#endif

namespace helicase {

// ─── Block ───────────────────────────────────────────────────────────────────
// A 64-byte aligned block returned by an input source.
// `ptr` is always a pointer to exactly 64 bytes (zero-padded for the final block).
// `size` is the number of valid bytes (1..64).
struct Block {
    const uint8_t* ptr = nullptr;
    size_t size = 0;
    explicit operator bool() const noexcept { return ptr != nullptr; }
};

// ─── SliceInput ──────────────────────────────────────────────────────────────
// Iterates a raw byte slice in 64-byte chunks.
// RANDOM_ACCESS = true: get_data() returns the full underlying slice.
class SliceInput {
    const uint8_t* data_;
    size_t size_;
    size_t pos_ = 0;
    alignas(64) uint8_t last_block_[64] = {};
    const uint8_t* cur_ptr_ = nullptr;
    size_t cur_size_ = 0;

public:
    static constexpr bool RANDOM_ACCESS = true;

    SliceInput(const uint8_t* data, size_t size) : data_(data), size_(size) {
        assert(size > 0);
        size_t tail = size % 64;
        size_t tail_start = (size / 64) * 64;
        if (tail > 0) std::memcpy(last_block_, data + tail_start, tail);
    }
    explicit SliceInput(const char* s) : SliceInput((const uint8_t*)s, std::strlen(s)) {}

    // Returns the next 64-byte block, or an empty Block at EOF.
    Block next() noexcept {
        size_t p = pos_;
        pos_ += 64;
        if (p + 64 <= size_) {
            cur_ptr_  = data_ + p;
            cur_size_ = 64;
        } else if (p < size_) {
            cur_ptr_  = last_block_;
            cur_size_ = size_ - p;
        } else {
            return {};
        }
        return {cur_ptr_, cur_size_};
    }

    // Returns the block most recently returned by next().
    const uint8_t* current_block() const noexcept { return cur_ptr_; }
    size_t current_block_size() const noexcept { return cur_size_; }

    // Full underlying data (only for RANDOM_ACCESS types).
    const uint8_t* get_data() const noexcept { return data_; }
    size_t get_size() const noexcept { return size_; }

    uint8_t first_byte() const noexcept { return data_[0]; }
};

// ─── MmapInput ───────────────────────────────────────────────────────────────
// Memory-mapped file; RANDOM_ACCESS = true.
#ifndef _WIN32
class MmapInput {
    const uint8_t* data_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;
    SliceInput* inner_ = nullptr;  // owned, heap-allocated

    alignas(64) uint8_t last_block_[64] = {};
    const uint8_t* cur_ptr_ = nullptr;
    size_t cur_size_ = 0;
    size_t pos_ = 0;

    void init() {
        size_t tail = size_ % 64;
        size_t tail_start = (size_ / 64) * 64;
        if (tail > 0) std::memcpy(last_block_, data_ + tail_start, tail);
    }

public:
    static constexpr bool RANDOM_ACCESS = true;

    explicit MmapInput(const std::string& path) {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) throw std::runtime_error("Cannot open file: " + path);
        struct stat st{};
        fstat(fd_, &st);
        size_ = (size_t)st.st_size;
        if (size_ == 0) { close(fd_); throw std::runtime_error("Empty file: " + path); }
        data_ = (const uint8_t*)mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (data_ == MAP_FAILED) { close(fd_); throw std::runtime_error("mmap failed: " + path); }
        init();
    }
    ~MmapInput() {
        if (data_) munmap((void*)data_, size_);
        if (fd_ >= 0) close(fd_);
    }
    MmapInput(const MmapInput&) = delete;
    MmapInput& operator=(const MmapInput&) = delete;
    MmapInput(MmapInput&& o) noexcept
        : data_(o.data_), size_(o.size_), fd_(o.fd_),
          cur_ptr_(o.cur_ptr_), cur_size_(o.cur_size_), pos_(o.pos_) {
        std::memcpy(last_block_, o.last_block_, 64);
        o.data_ = nullptr; o.fd_ = -1;
    }

    Block next() noexcept {
        size_t p = pos_;
        pos_ += 64;
        if (p + 64 <= size_) {
            cur_ptr_ = data_ + p; cur_size_ = 64;
        } else if (p < size_) {
            cur_ptr_ = last_block_; cur_size_ = size_ - p;
        } else {
            return {};
        }
        return {cur_ptr_, cur_size_};
    }

    const uint8_t* current_block() const noexcept { return cur_ptr_; }
    size_t current_block_size() const noexcept { return cur_size_; }
    const uint8_t* get_data() const noexcept { return data_; }
    size_t get_size() const noexcept { return size_; }
    uint8_t first_byte() const noexcept { return data_[0]; }
};
#endif // !_WIN32

// ─── FileInput ───────────────────────────────────────────────────────────────
// Buffered file reader (raw, no decompression). RANDOM_ACCESS = false.
// For compressed files, decompress externally before feeding.
class FileInput {
    static constexpr size_t BUF_SIZE = 1u << 16;  // 64 KB

    FILE* fp_ = nullptr;
    bool owns_ = false;
    std::vector<uint8_t> buf_;
    size_t buf_start_ = 0;  // offset of buf_[0] within the logical stream
    size_t buf_fill_ = 0;   // bytes currently in buf_
    size_t pos_ = 0;        // absolute position in stream
    bool eof_ = false;
    alignas(64) uint8_t block_[64] = {};
    const uint8_t* cur_ptr_ = nullptr;
    size_t cur_size_ = 0;
    uint8_t first_byte_ = 0;

    void refill() {
        // Shift unconsumed bytes to front
        size_t consumed = pos_ - buf_start_;
        if (consumed > 0 && consumed <= buf_fill_) {
            buf_fill_ -= consumed;
            std::memmove(buf_.data(), buf_.data() + consumed, buf_fill_);
            buf_start_ = pos_;
        }
        // Fill up
        size_t n = std::fread(buf_.data() + buf_fill_, 1, buf_.size() - buf_fill_, fp_);
        buf_fill_ += n;
        if (n == 0) eof_ = true;
    }

public:
    static constexpr bool RANDOM_ACCESS = false;

    explicit FileInput(const std::string& path) : buf_(BUF_SIZE) {
        fp_ = std::fopen(path.c_str(), "rb");
        if (!fp_) throw std::runtime_error("Cannot open file: " + path);
        owns_ = true;
        buf_fill_ = std::fread(buf_.data(), 1, buf_.size(), fp_);
        if (buf_fill_ == 0) throw std::runtime_error("Empty file: " + path);
        first_byte_ = buf_[0];
    }
    explicit FileInput(FILE* fp) : buf_(BUF_SIZE), fp_(fp), owns_(false) {
        buf_fill_ = std::fread(buf_.data(), 1, buf_.size(), fp_);
        if (buf_fill_ > 0) first_byte_ = buf_[0];
    }
    ~FileInput() { if (owns_ && fp_) std::fclose(fp_); }
    FileInput(const FileInput&) = delete;
    FileInput& operator=(const FileInput&) = delete;

    Block next() noexcept {
        // Need 64 bytes in buf starting at pos_
        size_t avail = buf_start_ + buf_fill_ - pos_;
        if (avail == 0) {
            if (eof_) return {};
            refill();
            avail = buf_start_ + buf_fill_ - pos_;
            if (avail == 0) return {};
        }
        if (avail >= 64) {
            cur_ptr_  = buf_.data() + (pos_ - buf_start_);
            cur_size_ = 64;
        } else {
            std::memcpy(block_, buf_.data() + (pos_ - buf_start_), avail);
            std::memset(block_ + avail, 0, 64 - avail);
            cur_ptr_  = block_;
            cur_size_ = avail;
            // Trigger EOF detection next call
        }
        size_t old_pos = pos_;
        pos_ += 64;
        // Refill if needed
        if (pos_ > buf_start_ + buf_fill_ && !eof_) refill();
        (void)old_pos;
        return {cur_ptr_, cur_size_};
    }

    const uint8_t* current_block() const noexcept { return cur_ptr_; }
    size_t current_block_size() const noexcept { return cur_size_; }
    uint8_t first_byte() const noexcept { return first_byte_; }

    // Not available for non-random-access inputs:
    // const uint8_t* get_data() -- not implemented
};

} // namespace helicase
