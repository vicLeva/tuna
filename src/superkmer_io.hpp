#pragma once

// On-disk / in-memory superkmer format:
//   [hdr_t len_bases][hdr_t min_pos][ceil(len/4) packed bytes]
//
// hdr_t is uint8_t when max superkmer length fits in 8 bits (2k-m ≤ 255),
// uint16_t otherwise.  The type is a compile-time constant deduced from k and m,
// so the header is as compact as possible without ever truncating a value.
//
//   Max superkmer length : 2k − m  (two k-mers sharing the boundary minimizer)
//   Max min_pos          : 2(k−m)  < 2k−m, so same type covers both fields
//
// For k ≤ 138 with m ≥ 21 (the common genomics range), hdr_t = uint8_t → 2-byte
// header per superkmer.  For larger k, hdr_t = uint16_t → 4-byte header.
//
// The sentinel value hdr_t::max() in min_pos means "no precomputed minimizer
// hash stored" — Phase 2 falls back to a full MinimizerWindow::reset().
//
// Packed encoding: base i is at bits 7-2*(i%4) of byte i/4, MSB-first.

#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <array>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <filesystem>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <unistd.h>

#ifdef TUNA_LZ4_BUCKETS
#  ifdef TUNA_HAVE_LZ4_H
#    include <lz4.h>
#  else
extern "C" {
int LZ4_compress_default(const char* src, char* dst, int src_size, int dst_capacity);
int LZ4_compress_fast(const char* src, char* dst, int src_size, int dst_capacity, int acceleration);
int LZ4_decompress_safe(const char* src, char* dst, int compressed_size, int dst_capacity);
int LZ4_compressBound(int input_size);
}
#  endif
#endif

// Returns the path to the superkmer file for partition p under work_dir.
inline std::string partition_path(const std::string& work_dir, size_t p)
{
    return work_dir + "hash_" + std::to_string(p) + ".superkmers";
}

class SuperkmerBucketFile {
public:
    SuperkmerBucketFile() = default;

    SuperkmerBucketFile(const SuperkmerBucketFile&) = delete;
    SuperkmerBucketFile& operator=(const SuperkmerBucketFile&) = delete;

    SuperkmerBucketFile(SuperkmerBucketFile&& o) noexcept
        : fd_(o.fd_), ok_(o.ok_)
    {
        o.fd_ = -1;
        o.ok_ = true;
    }

    SuperkmerBucketFile& operator=(SuperkmerBucketFile&& o) noexcept
    {
        if (this != &o) {
            close();
            fd_ = o.fd_;
            ok_ = o.ok_;
            o.fd_ = -1;
            o.ok_ = true;
        }
        return *this;
    }

    ~SuperkmerBucketFile() { close(); }

    void open(const std::string& path)
    {
        close();
#ifdef O_CLOEXEC
        constexpr int cloexec = O_CLOEXEC;
#else
        constexpr int cloexec = 0;
#endif
        fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | cloexec, 0666);
        ok_ = fd_ >= 0;
    }

    void open(const std::string& path, std::ios_base::openmode)
    {
        open(path);
    }

    void write(const char* data, std::streamsize n)
    {
        if (!ok_ || n < 0) {
            ok_ = false;
            return;
        }
        const char* cur = data;
        size_t left = static_cast<size_t>(n);
        while (left > 0) {
            const ssize_t written = ::write(fd_, cur, left);
            if (written <= 0) {
                ok_ = false;
                return;
            }
            cur += written;
            left -= static_cast<size_t>(written);
        }
    }

    void writev_pair(const char* a, size_t a_size, const char* b, size_t b_size)
    {
        if (!ok_) return;
        iovec iov[2] = {
            {const_cast<char*>(a), a_size},
            {const_cast<char*>(b), b_size}
        };
        size_t total = a_size + b_size;
        while (total > 0) {
            const ssize_t written = ::writev(fd_, iov, 2);
            if (written <= 0) {
                ok_ = false;
                return;
            }
            total -= static_cast<size_t>(written);
            size_t done = static_cast<size_t>(written);
            for (auto& v : iov) {
                if (done >= v.iov_len) {
                    done -= v.iov_len;
                    v.iov_base = static_cast<char*>(v.iov_base) + v.iov_len;
                    v.iov_len = 0;
                } else {
                    v.iov_base = static_cast<char*>(v.iov_base) + done;
                    v.iov_len -= done;
                    break;
                }
            }
        }
    }

    void close()
    {
        if (fd_ >= 0) {
            if (::close(fd_) != 0)
                ok_ = false;
            fd_ = -1;
        }
    }

    explicit operator bool() const noexcept { return ok_; }

private:
    int  fd_ = -1;
    bool ok_ = true;
};

namespace tuna_superkmer_detail {

#ifdef TUNA_LZ4_BUCKETS
inline constexpr char LZ4_BUCKET_MAGIC[] = {'T', 'U', 'N', 'A', 'L', 'Z', '4', 1};
inline constexpr size_t LZ4_BUCKET_MAGIC_SIZE = sizeof(LZ4_BUCKET_MAGIC);
inline constexpr uint8_t LZ4_BLOCK_RAW = 0;
inline constexpr uint8_t LZ4_BLOCK_COMPRESSED = 1;

inline void store_u32_le(char* b, uint32_t v)
{
    b[0] = static_cast<char>(v & 0xffu);
    b[1] = static_cast<char>((v >> 8) & 0xffu);
    b[2] = static_cast<char>((v >> 16) & 0xffu);
    b[3] = static_cast<char>((v >> 24) & 0xffu);
}

template <typename Out>
inline void write_block_payload(Out& out, const char* header, size_t header_size,
                                const char* payload, size_t payload_size)
{
    out.write(header, static_cast<std::streamsize>(header_size));
    out.write(payload, static_cast<std::streamsize>(payload_size));
}

inline void write_block_payload(SuperkmerBucketFile& out, const char* header, size_t header_size,
                                const char* payload, size_t payload_size)
{
    out.writev_pair(header, header_size, payload, payload_size);
}

inline uint32_t read_u32_le(const char* b) noexcept
{
    return static_cast<uint32_t>(static_cast<unsigned char>(b[0]))
         | (static_cast<uint32_t>(static_cast<unsigned char>(b[1])) << 8)
         | (static_cast<uint32_t>(static_cast<unsigned char>(b[2])) << 16)
         | (static_cast<uint32_t>(static_cast<unsigned char>(b[3])) << 24);
}

inline bool has_lz4_magic(const char* data) noexcept
{
    return std::memcmp(data, LZ4_BUCKET_MAGIC, LZ4_BUCKET_MAGIC_SIZE) == 0;
}

template <typename Out>
inline void write_lz4_bucket_magic(Out& out)
{
    out.write(LZ4_BUCKET_MAGIC, static_cast<std::streamsize>(LZ4_BUCKET_MAGIC_SIZE));
}
#endif

template <typename Out>
inline void write_superkmer_bucket_block(Out& out, const char* data, size_t size, bool compress = false)
{
    if (size == 0) return;
#ifdef TUNA_LZ4_BUCKETS
    if (!compress) {
        out.write(data, static_cast<std::streamsize>(size));
        return;
    }
    if (size > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("tuna: superkmer block too large for LZ4");

    const int raw_size = static_cast<int>(size);
    const int bound = LZ4_compressBound(raw_size);
    thread_local std::vector<char> compressed;
    compressed.resize(static_cast<size_t>(bound));
#ifdef TUNA_LZ4_ACCELERATION
    const int compressed_size = LZ4_compress_fast(
        data, compressed.data(), raw_size, bound, TUNA_LZ4_ACCELERATION);
#else
    const int compressed_size = LZ4_compress_default(data, compressed.data(), raw_size, bound);
#endif

    const bool use_compressed = compressed_size > 0
                             && static_cast<size_t>(compressed_size) + 9u < size;
    char header[9];
    store_u32_le(header, static_cast<uint32_t>(size));
    if (use_compressed) {
        store_u32_le(header + 4, static_cast<uint32_t>(compressed_size));
        header[8] = static_cast<char>(LZ4_BLOCK_COMPRESSED);
        write_block_payload(out, header, sizeof(header), compressed.data(),
                            static_cast<size_t>(compressed_size));
    } else {
        store_u32_le(header + 4, static_cast<uint32_t>(size));
        header[8] = static_cast<char>(LZ4_BLOCK_RAW);
        write_block_payload(out, header, sizeof(header), data, size);
    }
#else
    out.write(data, static_cast<std::streamsize>(size));
#endif
}

#ifdef TUNA_LZ4_BUCKETS
inline bool load_lz4_bucket(const std::string& path, std::vector<char>& out_data)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;

    char magic[LZ4_BUCKET_MAGIC_SIZE];
    in.read(magic, static_cast<std::streamsize>(LZ4_BUCKET_MAGIC_SIZE));
    if (in.gcount() != static_cast<std::streamsize>(LZ4_BUCKET_MAGIC_SIZE))
        return false;
    if (!has_lz4_magic(magic))
        return false;

    std::error_code ec;
    const auto compressed_size = std::filesystem::file_size(path, ec);
    if (!ec)
        out_data.reserve(static_cast<size_t>(compressed_size) * 2u);

    char header[9];
    std::vector<char> stored;
    while (true) {
        in.read(header, sizeof(header));
        if (in.gcount() == 0 && in.eof()) break;
        if (in.gcount() != static_cast<std::streamsize>(sizeof(header)))
            throw std::runtime_error("tuna: truncated LZ4 superkmer bucket: " + path);

        const uint32_t raw_size = read_u32_le(header);
        const uint32_t stored_size = read_u32_le(header + 4);
        const uint8_t flags = static_cast<uint8_t>(header[8]);
        if (raw_size == 0 || stored_size == 0)
            throw std::runtime_error("tuna: invalid LZ4 superkmer bucket block: " + path);

        const size_t old_size = out_data.size();
        out_data.resize(old_size + raw_size);
        if (flags == LZ4_BLOCK_RAW) {
            if (stored_size != raw_size)
                throw std::runtime_error("tuna: invalid raw LZ4 superkmer bucket block: " + path);
            in.read(out_data.data() + old_size, static_cast<std::streamsize>(raw_size));
            if (!in)
                throw std::runtime_error("tuna: truncated raw LZ4 superkmer bucket block: " + path);
        } else if (flags == LZ4_BLOCK_COMPRESSED) {
            stored.resize(stored_size);
            in.read(stored.data(), static_cast<std::streamsize>(stored_size));
            if (!in)
                throw std::runtime_error("tuna: truncated compressed LZ4 superkmer bucket block: " + path);
            const int decoded = LZ4_decompress_safe(
                stored.data(), out_data.data() + old_size,
                static_cast<int>(stored_size), static_cast<int>(raw_size));
            if (decoded != static_cast<int>(raw_size))
                throw std::runtime_error("tuna: failed to decode LZ4 superkmer bucket block: " + path);
        } else {
            throw std::runtime_error("tuna: unknown LZ4 superkmer bucket block type: " + path);
        }
    }

    return true;
}
#endif

} // namespace tuna_superkmer_detail

// Deduce the header integer type for a given (k, m) pair.
// uint8_t  when 2k − m ≤ 255  (header fits in one byte, 2-byte header total)
// uint16_t otherwise           (4-byte header, required for k ≥ 139 at m=21)
template <uint16_t k, uint16_t m>
using sk_hdr_t = std::conditional_t<(2u * k - m <= 255u), uint8_t, uint16_t>;

// Sentinel stored in min_pos when no precomputed minimizer hash is available.
template <uint16_t k, uint16_t m>
static constexpr sk_hdr_t<k, m> sk_no_min = std::numeric_limits<sk_hdr_t<k, m>>::max();


// ─── Writer ───────────────────────────────────────────────────────────────────
//
// SuperkmerWriter uses a raw POD buffer instead of std::string to avoid
// zero-initialisation on resize.  std::string::resize zeroes new bytes even when
// immediately overwritten; this buffer grows by doubling but skips zero-fill.

struct SuperkmerWriteBlock
{
    size_t partition = 0;
    char*  data      = nullptr;
    size_t size      = 0;
    size_t capacity  = 0;

    SuperkmerWriteBlock() = default;
    SuperkmerWriteBlock(size_t p, char* d, size_t s, size_t c) noexcept
        : partition(p), data(d), size(s), capacity(c) {}

    SuperkmerWriteBlock(const SuperkmerWriteBlock&) = delete;
    SuperkmerWriteBlock& operator=(const SuperkmerWriteBlock&) = delete;

    SuperkmerWriteBlock(SuperkmerWriteBlock&& o) noexcept
        : partition(o.partition), data(o.data), size(o.size), capacity(o.capacity)
    { o.data = nullptr; o.size = 0; o.capacity = 0; }

    SuperkmerWriteBlock& operator=(SuperkmerWriteBlock&& o) noexcept
    {
        if (this != &o) {
            ::operator delete(data);
            partition = o.partition;
            data      = o.data;
            size      = o.size;
            capacity  = o.capacity;
            o.data = nullptr;
            o.size = 0;
            o.capacity = 0;
        }
        return *this;
    }

    ~SuperkmerWriteBlock() { ::operator delete(data); }
};

template <uint16_t k, uint16_t m>
struct SuperkmerWriter
{
    using hdr_t = sk_hdr_t<k, m>;              // superkmer header type (local alias)
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);  // header bytes per record

    char*  raw_  = nullptr;  // raw buffer pointer
    size_t sz_   = 0;        // used bytes
    size_t cap_  = 0;        // allocated bytes
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
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst,                &len,     sizeof(hdr_t));
        std::memcpy(dst + sizeof(hdr_t), &min_pos, sizeof(hdr_t));
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
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst,                 &len,     sizeof(hdr_t));
        std::memcpy(dst + sizeof(hdr_t), &min_pos, sizeof(hdr_t));
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

    template <typename GetNtBase>
    void append_nt(GetNtBase&& get_nt, size_t start, hdr_t len, hdr_t min_pos)
    {
        const size_t packed_bytes = (len + 3u) / 4u;
        char* dst = reserve_inline(HDR_BYTES + packed_bytes);
        std::memcpy(dst,                 &len,     sizeof(hdr_t));
        std::memcpy(dst + sizeof(hdr_t), &min_pos, sizeof(hdr_t));
        uint8_t* packed = reinterpret_cast<uint8_t*>(dst + HDR_BYTES);

        size_t i = 0;
        for (; i + 4 <= static_cast<size_t>(len); i += 4) {
            const uint8_t nt0 = get_nt(start + i);
            const uint8_t nt1 = get_nt(start + i + 1);
            const uint8_t nt2 = get_nt(start + i + 2);
            const uint8_t nt3 = get_nt(start + i + 3);
            const uint8_t b0 = nt0 ^ (nt0 >> 1);
            const uint8_t b1 = nt1 ^ (nt1 >> 1);
            const uint8_t b2 = nt2 ^ (nt2 >> 1);
            const uint8_t b3 = nt3 ^ (nt3 >> 1);
            packed[i >> 2] = static_cast<uint8_t>((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
        }
        if (i < static_cast<size_t>(len)) {
            uint8_t tail = 0;
            for (size_t j = i; j < static_cast<size_t>(len); ++j) {
                const uint8_t nt = get_nt(start + j);
                const uint8_t b = nt ^ (nt >> 1);
                tail |= static_cast<uint8_t>(b << (6u - 2u * (j - i)));
            }
            packed[i >> 2] = tail;
        }
    }

    bool needs_flush() const noexcept { return sz_ >= flush_threshold; }

    SuperkmerWriteBlock release_block(size_t partition)
    {
        SuperkmerWriteBlock block(partition, raw_, sz_, cap_);
        cap_ = flush_threshold;
        raw_ = static_cast<char*>(::operator new(cap_));
        sz_  = 0;
        return block;
    }

    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (sz_ == 0) return;
        std::lock_guard<std::mutex> g(mtx);
        tuna_superkmer_detail::write_superkmer_bucket_block(file, raw_, sz_);
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
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);
    static constexpr size_t MAX_LEN = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr size_t MAX_RECORD_BYTES = HDR_BYTES + (MAX_LEN + 3u) / 4u;

    explicit SuperkmerReader(const std::string& path)
    {
        std::error_code ec;
        const auto fsz = std::filesystem::file_size(path, ec);
        if (ec || fsz == 0) return;
#ifdef TUNA_LZ4_BUCKETS
        if (tuna_superkmer_detail::load_lz4_bucket(path, data_)) {
            size_ = data_.size();
            cur_ = data_.data();
            end_ = data_.data() + data_.size();
            ok_ = true;
            return;
        }
#endif
        size_ = static_cast<size_t>(fsz);
        data_.resize(size_);
        std::ifstream in(path, std::ios::binary);
        if (!in) return;
        in.read(data_.data(), static_cast<std::streamsize>(data_.size()));
        if (!in) {
            data_.clear();
            size_ = 0;
            return;
        }
        cur_ = data_.data();
        end_ = data_.data() + data_.size();
        ok_ = true;
    }

    ~SuperkmerReader() = default;

    SuperkmerReader(const SuperkmerReader&)            = delete;
    SuperkmerReader& operator=(const SuperkmerReader&) = delete;

    bool next()
    {
        if (!ok_) return false;
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        hdr_t len, mp;
        std::memcpy(&len, cur_,                 sizeof(hdr_t));
        std::memcpy(&mp,  cur_ + sizeof(hdr_t), sizeof(hdr_t));
        if (len == 0) return false;
        min_pos_ = mp;

        const size_t packed_bytes = (static_cast<size_t>(len) + 3u) / 4u;
        if (HDR_BYTES + packed_bytes > MAX_RECORD_BYTES) return false;
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES + packed_bytes) > end_) return false;

        record_      = cur_;
        record_size_ = HDR_BYTES + packed_bytes;
        ptr_         = reinterpret_cast<const uint8_t*>(cur_ + HDR_BYTES);
        len_         = len;
        cur_        += record_size_;
        return true;
    }

    const uint8_t* packed_data() const { return ptr_; }
    size_t         size()        const { return len_; }
    hdr_t          min_pos()     const { return min_pos_; }
    size_t         file_size()   const { return size_; }
    bool           ok()          const { return ok_; }
    const char*    record_data() const { return record_; }
    size_t         record_size() const { return record_size_; }

private:
    std::vector<char> data_;
    bool              ok_          = false;
    size_t            size_        = 0;       // file size in bytes
    const char*       cur_         = nullptr;
    const char*       end_         = nullptr;
    const char*       record_      = nullptr;
    size_t            record_size_ = 0;
    const uint8_t*    ptr_         = nullptr; // packed data of current superkmer
    size_t            len_         = 0;       // length of current superkmer (bases)
    hdr_t             min_pos_     = 0;       // minimizer position in current superkmer
};


// ─── In-memory reader (streaming mode) ───────────────────────────────────────

template <uint16_t k, uint16_t m>
struct MemoryReader
{
    using hdr_t = sk_hdr_t<k, m>;
    static constexpr size_t HDR_BYTES = 2 * sizeof(hdr_t);

    MemoryReader() = default;
    explicit MemoryReader(const std::string& data) noexcept
        : cur_(data.data()), end_(data.data() + data.size()) {}

    bool next() noexcept
    {
        if (cur_ + static_cast<ptrdiff_t>(HDR_BYTES) > end_) return false;
        hdr_t len, mp;
        std::memcpy(&len, cur_,                sizeof(hdr_t));
        std::memcpy(&mp,  cur_ + sizeof(hdr_t), sizeof(hdr_t));
        if (len == 0) return false;
        min_pos_ = mp;
        cur_ += HDR_BYTES;
        const size_t packed_bytes = (static_cast<size_t>(len) + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(packed_bytes) > end_) return false;
        ptr_  = reinterpret_cast<const uint8_t*>(cur_);
        len_  = len;
        cur_ += packed_bytes;
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
