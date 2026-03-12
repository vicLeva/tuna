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
// Dispatch:
//   Plain files  — zero-copy mmap + helicase SIMD (split_non_actg).
//   GZ files     — 64 MB streaming window (SeqReader) + split_actg() on each sequence.
//
// The internal gz buffer only grows, so reusing one SeqSource across many files
// eliminates per-file allocation overhead.

#include "nt_hash.hpp"
#include <helicase.hpp>
#include <zlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <cctype>
#include <algorithm>


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


// ── SeqReader (internal) ──────────────────────────────────────────────────────
// Streaming FASTA/FASTQ reader with optional gzip decompression.
// Used by SeqSource for gz files; exposed publicly for callers that need the
// sequence-at-a-time interface (e.g. partition_kmc pre-scan helper).

enum class SeqFormat { FASTA, FASTQ };

struct SeqReader
{
    SeqReader()
    {
        buf_.reserve(6 << 20);
        seq_.reserve(6 << 20);
    }

    ~SeqReader() { close(); }

    bool load(const std::string& path)
    {
        close();
        std::string inner = path;
        bool gz = false;
        if (inner.size() > 3) {
            auto tail = inner.substr(inner.size() - 3);
            for (char& c : tail) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (tail == ".gz") { gz = true; inner.resize(inner.size() - 3); }
        }
        fmt_ = detect_format(inner);
        return gz ? load_gz(path) : load_plain(path);
    }

    bool read_next_seq()
    {
        if (gz_ && !gz_eof_ && fmt_ == SeqFormat::FASTQ) {
            if (buf_.size() - rpos_ < (4u << 20))
                gz_slide_and_refill();
        }
        seq_.clear();
        return fmt_ == SeqFormat::FASTQ ? read_fastq() : read_fasta();
    }

    const char* seq()     const { return seq_.data(); }
    size_t      seq_len() const { return seq_.size(); }

    void close()
    {
        if (gz_) { gzclose(gz_); gz_ = nullptr; gz_eof_ = false; }
    }

private:
    SeqFormat         fmt_ = SeqFormat::FASTA;
    std::vector<char> buf_;
    std::vector<char> seq_;
    size_t            rpos_ = 0;
    gzFile gz_     = nullptr;
    bool   gz_eof_ = false;
    static constexpr size_t GZ_CHUNK = 64u << 20;

    bool load_plain(const std::string& path)
    {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) return false;
        const auto sz = static_cast<size_t>(f.tellg());
        if (sz == 0) return false;
        f.seekg(0);
        buf_.resize(sz);
        f.read(buf_.data(), static_cast<std::streamsize>(sz));
        rpos_ = 0;
        return true;
    }

    bool load_gz(const std::string& path)
    {
        gz_ = gzopen(path.c_str(), "rb");
        if (!gz_) return false;
        gzbuffer(gz_, 256u << 10);
        gz_eof_ = false;
        buf_.clear();
        rpos_ = 0;
        buf_.resize(GZ_CHUNK);
        const size_t got = gz_read_into(0, GZ_CHUNK);
        buf_.resize(got);
        return got > 0;
    }

    size_t gz_read_into(size_t offset, size_t ask)
    {
        size_t got = 0;
        while (got < ask && !gz_eof_) {
            const unsigned chunk = static_cast<unsigned>(std::min(ask - got, size_t(1u << 30)));
            const int n = gzread(gz_, buf_.data() + offset + got, chunk);
            if (n <= 0) { gz_eof_ = true; break; }
            got += static_cast<size_t>(n);
        }
        return got;
    }

    void gz_slide_and_refill()
    {
        if (!gz_ || gz_eof_) return;
        const size_t rem = (rpos_ < buf_.size()) ? buf_.size() - rpos_ : 0;
        if (rem > 0 && rpos_ > 0)
            std::memmove(buf_.data(), buf_.data() + rpos_, rem);
        rpos_ = 0;
        buf_.resize(rem + GZ_CHUNK);
        const size_t got = gz_read_into(rem, GZ_CHUNK);
        buf_.resize(rem + got);
    }

    bool read_fasta()
    {
        while (rpos_ < buf_.size() && buf_[rpos_] != '>') ++rpos_;
        if (rpos_ >= buf_.size()) {
            if (!gz_ || gz_eof_) return false;
            gz_slide_and_refill();
            while (rpos_ < buf_.size() && buf_[rpos_] != '>') ++rpos_;
            if (rpos_ >= buf_.size()) return false;
        }
        ++rpos_;

        while (true) {
            const char* nl = static_cast<const char*>(
                std::memchr(buf_.data() + rpos_, '\n', buf_.size() - rpos_));
            if (nl) { rpos_ = static_cast<size_t>(nl - buf_.data()) + 1; break; }
            if (!gz_ || gz_eof_) return false;
            gz_slide_and_refill();
            if (rpos_ >= buf_.size()) return false;
        }

        while (true) {
            if (rpos_ >= buf_.size()) {
                if (!gz_ || gz_eof_) break;
                gz_slide_and_refill();
                if (rpos_ >= buf_.size()) break;
            }
            if (buf_[rpos_] == '>') break;

            const char*  start = buf_.data() + rpos_;
            const size_t rem   = buf_.size() - rpos_;
            const char*  nl    = static_cast<const char*>(std::memchr(start, '\n', rem));

            if (!nl) {
                if (gz_ && !gz_eof_) {
                    size_t line_len = rem;
                    if (line_len > 0 && start[line_len - 1] == '\r') --line_len;
                    if (line_len > 0) {
                        const size_t old_sz = seq_.size();
                        seq_.resize(old_sz + line_len);
                        std::memcpy(seq_.data() + old_sz, start, line_len);
                    }
                    rpos_ = buf_.size();
                    gz_slide_and_refill();
                    continue;
                }
                size_t line_len = rem;
                if (line_len > 0 && start[line_len - 1] == '\r') --line_len;
                if (line_len > 0) {
                    const size_t old_sz = seq_.size();
                    seq_.resize(old_sz + line_len);
                    std::memcpy(seq_.data() + old_sz, start, line_len);
                }
                rpos_ = buf_.size();
                break;
            }

            size_t line_len = static_cast<size_t>(nl - start);
            if (line_len > 0 && start[line_len - 1] == '\r') --line_len;
            if (line_len > 0) {
                const size_t old_sz = seq_.size();
                seq_.resize(old_sz + line_len);
                std::memcpy(seq_.data() + old_sz, start, line_len);
            }
            rpos_ += static_cast<size_t>(nl - start) + 1;
        }

        return !seq_.empty();
    }

    bool read_fastq()
    {
        if (rpos_ >= buf_.size()) return false;
        const char* end = buf_.data() + buf_.size();
        const char* p   = buf_.data() + rpos_;
        const char* nl;

        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        if (!nl) return false;
        p = nl + 1;

        const char* seq_start = p;
        nl = static_cast<const char*>(std::memchr(seq_start, '\n', end - seq_start));
        size_t slen = nl ? static_cast<size_t>(nl - seq_start) : static_cast<size_t>(end - seq_start);
        if (slen > 0 && seq_start[slen - 1] == '\r') --slen;
        if (slen > 0) seq_.assign(seq_start, seq_start + slen);
        p = nl ? nl + 1 : end;

        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        p = nl ? nl + 1 : end;

        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        rpos_ = nl ? static_cast<size_t>(nl - buf_.data()) + 1 : buf_.size();

        return !seq_.empty();
    }

    static SeqFormat detect_format(const std::string& name)
    {
        auto has_suffix = [&](const char* suf) -> bool {
            const size_t n = std::strlen(suf);
            if (name.size() < n) return false;
            return std::equal(name.end() - static_cast<std::ptrdiff_t>(n), name.end(),
                              suf, [](char a, char b) {
                                  return std::tolower(static_cast<unsigned char>(a)) == b;
                              });
        };
        if (has_suffix(".fastq") || has_suffix(".fq")) return SeqFormat::FASTQ;
        return SeqFormat::FASTA;
    }
};


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
        if (is_gz) process_gz(path, std::forward<F>(on_chunk));
        else       process_plain(path, std::forward<F>(on_chunk));
    }

private:
    SeqReader gz_reader_; // reused across gz files; buffer only grows

    template <typename F>
    void process_gz(const std::string& path, F&& on_chunk)
    {
        if (!gz_reader_.load(path)) return;
        while (gz_reader_.read_next_seq())
            split_actg(gz_reader_.seq(), gz_reader_.seq_len(),
                       std::forward<F>(on_chunk));
    }

    template <typename F>
    void process_plain(const std::string& path, F&& on_chunk)
    {
        helicase::MmapInput inp(path);
        if (inp.first_byte() == '@') {
            helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
            while (p.next()) {
                auto [ptr, len] = p.get_dna_raw();
                on_chunk(ptr, len);
            }
        } else {
            helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
            while (p.next()) {
                auto [ptr, len] = p.get_dna_raw();
                on_chunk(ptr, len);
            }
        }
    }

    // Split a sequence on non-ACTG characters and call on_chunk for each run.
    template <typename F>
    static void split_actg(const char* seq, size_t len, F&& on_chunk)
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
};
