#pragma once

// SeqSource — unified sequence input over plain and gzip FASTA/FASTQ.
//
// Single interface for all Phase 1 readers:
//
//   SeqSource src;                          // one instance per thread, reuse across files
//   src.process(path, [](const char* chunk, size_t len) {
//       // chunk is ACTG-only, newline-free; len >= 1
//   });
//
// All input types (plain/gz, FASTA/FASTQ) are handled uniformly by helicase SIMD
// parsers.  Plain files use zero-copy mmap (MmapInput).  GZ files stream through
// GzInput, a zlib-backed Block-producer that satisfies the helicase input interface.
//
// GzInput is also usable directly — e.g. in producer-consumer harnesses that need
// to decompress and parse independently from SeqSource.

#include "nt_hash.hpp"
#include <helicase.hpp>
#include <zlib.h>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <stdexcept>


// ── helicase config ───────────────────────────────────────────────────────────
// DNA string, ignore headers, split at every non-ACTG character.
// Each p.next() / p.get_dna_raw() delivers one ACTG-only, newline-free run.

static constexpr helicase::Config HELICASE_ACTG =
    helicase::ParserOptions()
        .ignore_headers()
        .dna_string()
        .split_non_actg()
        .config()
    & ~helicase::advanced::RETURN_RECORD;


// ── GzInput ───────────────────────────────────────────────────────────────────
// helicase-compatible Block-producer backed by zlib gz decompression.
// RANDOM_ACCESS = false; use with FastaParser<CFG, GzInput> / FastqParser<CFG, GzInput>.
// The constructor eagerly reads the first block so first_byte() is available before
// any call to next().

class GzInput {
    static constexpr size_t BUF_SIZE = 1u << 16;   // 64 KB ring buffer

    gzFile gz_   = nullptr;
    bool   eof_  = false;
    std::vector<uint8_t> buf_;
    size_t buf_start_ = 0;
    size_t buf_fill_  = 0;
    size_t pos_       = 0;
    alignas(64) uint8_t block_[64] = {};
    const uint8_t* cur_ptr_    = nullptr;
    size_t         cur_size_   = 0;
    uint8_t        first_byte_ = 0;

    void refill() {
        size_t consumed = pos_ - buf_start_;
        if (consumed > 0 && consumed <= buf_fill_) {
            buf_fill_ -= consumed;
            std::memmove(buf_.data(), buf_.data() + consumed, buf_fill_);
            buf_start_ = pos_;
        }
        while (buf_fill_ < buf_.size() && !eof_) {
            int n = gzread(gz_, buf_.data() + buf_fill_,
                           static_cast<unsigned>(buf_.size() - buf_fill_));
            if (n <= 0) { eof_ = true; break; }
            buf_fill_ += static_cast<size_t>(n);
        }
    }

public:
    static constexpr bool RANDOM_ACCESS = false;

    explicit GzInput(const std::string& path) : buf_(BUF_SIZE) {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) throw std::runtime_error("Cannot open: " + path);
        gzbuffer(gz_, 256u << 10);
        int n = gzread(gz_, buf_.data(), static_cast<unsigned>(buf_.size()));
        if (n <= 0) { gzclose(gz_); gz_ = nullptr; throw std::runtime_error("Cannot read: " + path); }
        buf_fill_   = static_cast<size_t>(n);
        first_byte_ = buf_[0];
        // eof_ is detected lazily by refill() when gzread returns 0.
    }
    ~GzInput() { if (gz_) gzclose(gz_); }
    GzInput(const GzInput&) = delete;
    GzInput& operator=(const GzInput&) = delete;
    GzInput(GzInput&& o) noexcept
        : gz_(o.gz_), eof_(o.eof_), buf_(std::move(o.buf_)),
          buf_start_(o.buf_start_), buf_fill_(o.buf_fill_), pos_(o.pos_),
          cur_ptr_(o.cur_ptr_), cur_size_(o.cur_size_), first_byte_(o.first_byte_) {
        std::memcpy(block_, o.block_, 64);
        o.gz_ = nullptr;
    }

    helicase::Block next() noexcept {
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
            pos_ += 64;
            // Lookahead: prefill the buffer before the next next() call.
            if (pos_ > buf_start_ + buf_fill_ && !eof_) refill();
        } else {
            // Last (partial) block.  Advance pos_ to the exact end of data so
            // that the next next() call gets avail==0 and returns {} correctly,
            // avoiding unsigned underflow if pos_ were advanced by a full 64.
            std::memcpy(block_, buf_.data() + (pos_ - buf_start_), avail);
            std::memset(block_ + avail, 0, 64 - avail);
            cur_ptr_  = block_;
            cur_size_ = avail;
            pos_ = buf_start_ + buf_fill_;
        }
        return {cur_ptr_, cur_size_};
    }

    const uint8_t* current_block()      const noexcept { return cur_ptr_; }
    size_t         current_block_size() const noexcept { return cur_size_; }
    uint8_t        first_byte()         const noexcept { return first_byte_; }
};


// ── split_actg (free function) ────────────────────────────────────────────────
// Split a raw sequence on non-ACTG characters and call on_chunk(ptr, len) for
// each ACTG-only run.  Header-free, usable outside SeqSource.

template <typename F>
inline void split_actg(const char* seq, size_t len, F&& on_chunk)
{
    size_t start  = 0;
    bool   in_run = false;
    for (size_t i = 0; i < len; ++i) {
        if (nt_hash::is_dna(seq[i])) {
            if (!in_run) { start = i; in_run = true; }
        } else {
            if (in_run) { on_chunk(seq + start, i - start); in_run = false; }
        }
    }
    if (in_run) on_chunk(seq + start, len - start);
}


// ── Dispatch helper ───────────────────────────────────────────────────────────
// Detect FASTA vs FASTQ from the first byte, instantiate the correct helicase
// parser, and call on_chunk for every ACTG-only run.  Works for any helicase
// input type I (MmapInput, GzInput, …).

template <typename I, typename F>
inline void helicase_parse(I inp, F&& on_chunk)
{
    if (inp.first_byte() == '@') {
        helicase::FastqParser<HELICASE_ACTG, I> p(std::move(inp));
        while (p.next()) {
            auto [ptr, len] = p.get_dna_raw();
            on_chunk(ptr, len);
        }
    } else {
        helicase::FastaParser<HELICASE_ACTG, I> p(std::move(inp));
        while (p.next()) {
            auto [ptr, len] = p.get_dna_raw();
            on_chunk(ptr, len);
        }
    }
}


// ── SeqSource ─────────────────────────────────────────────────────────────────

struct SeqSource
{
    SeqSource() = default;

    // Call on_chunk(const char* chunk, size_t len) for every ACTG-only run in path.
    template <typename F>
    void process(const std::string& path, F&& on_chunk)
    {
        const bool is_gz = path.size() > 3 &&
                           path.compare(path.size() - 3, 3, ".gz") == 0;
        if (is_gz)
            helicase_parse(GzInput(path), std::forward<F>(on_chunk));
        else
            helicase_parse(helicase::MmapInput(path), std::forward<F>(on_chunk));
    }
};
