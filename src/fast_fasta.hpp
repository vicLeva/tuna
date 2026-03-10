#pragma once

// SeqReader — unified FASTA/FASTQ reader with optional gzip decompression.
//
// Drop-in replacement for the old FastFastaReader (now supports FASTQ and .gz).
//
// Design: one instance per worker thread, kept alive across files.  Call load()
// for each new file; the internal buffers only grow, eliminating per-file
// allocation overhead.  Sequences are returned via seq() / seq_len() and are
// valid until the next call to read_next_seq() or load().
//
// Format detection (from file extension, case-insensitive; .gz stripped first):
//   FASTA : .fa  .fna  .fasta  (and any other extension)
//   FASTQ : .fq  .fastq
//
// Gzip decompression: streaming chunked reads — at most GZ_CHUNK bytes are
// held in buf_ at any time.  This avoids OOM on large compressed files (e.g.
// 30 GB decompressed FASTQ).  The sliding-window refill is transparent to the
// FASTA/FASTQ parsers via gz_slide_and_refill().

#include <zlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <cctype>
#include <algorithm>


enum class SeqFormat { FASTA, FASTQ };

struct SeqReader
{
    SeqReader()
    {
        // Pre-allocate for a typical genome file (~5 MB); grows as needed.
        buf_.reserve(6 << 20);
        seq_.reserve(6 << 20);
    }

    ~SeqReader() { close(); }

    // Load a new file.  Auto-detects format and compression from extension.
    // Reuses existing buffer capacity.  Returns false if the file cannot be
    // opened or yields no data.
    bool load(const std::string& path)
    {
        close(); // close any previously open gz handle

        // Strip .gz suffix to reveal inner extension.
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

    // Advance to the next sequence in the current file.
    // Returns false when no more sequences are available.
    bool read_next_seq()
    {
        // For FASTQ streaming: proactively refill before parsing each record
        // so that the read_fastq() parser always sees a complete record.
        if (gz_ && !gz_eof_ && fmt_ == SeqFormat::FASTQ) {
            if (buf_.size() - rpos_ < (4u << 20)) // < 4 MB remaining
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
    std::vector<char> buf_;      // raw (or decompressed window) file content
    std::vector<char> seq_;      // current stripped sequence; reused
    size_t            rpos_ = 0; // read cursor in buf_

    // Streaming gz state
    gzFile gz_     = nullptr;
    bool   gz_eof_ = false;
    static constexpr size_t GZ_CHUNK = 64u << 20; // 64 MB window


    // ── Loaders ───────────────────────────────────────────────────────────────

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
        gzbuffer(gz_, 256u << 10); // 256 KB zlib internal buffer
        gz_eof_ = false;
        buf_.clear();
        rpos_ = 0;

        // Fill first chunk
        buf_.resize(GZ_CHUNK);
        const size_t got = gz_read_into(0, GZ_CHUNK);
        buf_.resize(got);
        return got > 0;
    }

    // Read up to `ask` bytes from gz_ into buf_[offset..].
    // buf_ must already be sized to at least offset + ask.
    // Returns bytes actually read; sets gz_eof_ on stream end.
    size_t gz_read_into(size_t offset, size_t ask)
    {
        size_t got = 0;
        while (got < ask && !gz_eof_) {
            const unsigned chunk = static_cast<unsigned>(
                std::min(ask - got, size_t(1u << 30)));
            const int n = gzread(gz_, buf_.data() + offset + got, chunk);
            if (n <= 0) { gz_eof_ = true; break; }
            got += static_cast<size_t>(n);
        }
        return got;
    }

    // Slide [rpos_, buf_.size()) to the front of buf_ and append up to
    // GZ_CHUNK more bytes from the gz stream.  Called when the parser is
    // running low on data.
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


    // ── Parsers ───────────────────────────────────────────────────────────────

    // FASTA: collect sequence lines between '>' headers, stripping newlines.
    // Handles streaming gz files: refills the window when the buffer is
    // exhausted mid-record.
    bool read_fasta()
    {
        // Seek to the next '>'
        while (rpos_ < buf_.size() && buf_[rpos_] != '>') ++rpos_;
        if (rpos_ >= buf_.size()) {
            if (!gz_ || gz_eof_) return false;
            gz_slide_and_refill();
            while (rpos_ < buf_.size() && buf_[rpos_] != '>') ++rpos_;
            if (rpos_ >= buf_.size()) return false;
        }
        ++rpos_; // consume '>'

        // Skip header line
        while (true) {
            const char* nl = static_cast<const char*>(
                std::memchr(buf_.data() + rpos_, '\n', buf_.size() - rpos_));
            if (nl) {
                rpos_ = static_cast<size_t>(nl - buf_.data()) + 1;
                break;
            }
            if (!gz_ || gz_eof_) return false; // malformed: no newline in header
            gz_slide_and_refill();
            if (rpos_ >= buf_.size()) return false;
        }

        // Accumulate sequence lines until the next '>' or EOF.
        while (true) {
            // Refill if we have exhausted the current window.
            if (rpos_ >= buf_.size()) {
                if (!gz_ || gz_eof_) break; // EOF — last sequence
                gz_slide_and_refill();
                if (rpos_ >= buf_.size()) break; // still empty → EOF
            }

            if (buf_[rpos_] == '>') break; // next sequence starts here

            const char*  start = buf_.data() + rpos_;
            const size_t rem   = buf_.size() - rpos_;
            const char*  nl    = static_cast<const char*>(std::memchr(start, '\n', rem));

            if (!nl) {
                // No newline visible in the current window.
                if (gz_ && !gz_eof_) {
                    // Append the partial line fragment then slide-refill.
                    size_t line_len = rem;
                    if (line_len > 0 && start[line_len - 1] == '\r') --line_len;
                    if (line_len > 0) {
                        const size_t old_sz = seq_.size();
                        seq_.resize(old_sz + line_len);
                        std::memcpy(seq_.data() + old_sz, start, line_len);
                    }
                    rpos_ = buf_.size(); // mark as consumed before slide
                    gz_slide_and_refill();
                    continue;
                }
                // Plain file or gz at EOF — last line without newline.
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

            // Normal case: found a newline within the current window.
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

    // FASTQ: 4-line records — @header / sequence / +ignored / quality.
    // Only the sequence line is kept; quality is skipped with a single memchr.
    // For gz files, gz_slide_and_refill() is called from read_next_seq() before
    // each record, so a complete record is always available in buf_.
    bool read_fastq()
    {
        if (rpos_ >= buf_.size()) return false;
        const char* end = buf_.data() + buf_.size();
        const char* p   = buf_.data() + rpos_;
        const char* nl;

        // Line 1: @header — skip entirely.
        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        if (!nl) return false;
        p = nl + 1;

        // Line 2: sequence.
        const char* seq_start = p;
        nl = static_cast<const char*>(std::memchr(seq_start, '\n', end - seq_start));
        size_t slen = nl ? static_cast<size_t>(nl - seq_start)
                         : static_cast<size_t>(end - seq_start);
        if (slen > 0 && seq_start[slen - 1] == '\r') --slen;
        if (slen > 0)
            seq_.assign(seq_start, seq_start + slen);
        p = nl ? nl + 1 : end;

        // Line 3: '+' separator — skip.
        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        p = nl ? nl + 1 : end;

        // Line 4: quality string — skip, advance cursor past it.
        nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
        rpos_ = nl ? static_cast<size_t>(nl - buf_.data()) + 1 : buf_.size();

        return !seq_.empty();
    }


    // ── Extension helpers ─────────────────────────────────────────────────────

    static SeqFormat detect_format(const std::string& name)
    {
        // Check whether `name` ends with `suf` (case-insensitive).
        auto has_suffix = [&](const char* suf) -> bool {
            const size_t n = std::strlen(suf);
            if (name.size() < n) return false;
            return std::equal(name.end() - static_cast<std::ptrdiff_t>(n), name.end(),
                              suf, [](char a, char b) {
                                  return std::tolower(static_cast<unsigned char>(a)) == b;
                              });
        };
        if (has_suffix(".fastq") || has_suffix(".fq")) return SeqFormat::FASTQ;
        return SeqFormat::FASTA; // .fa / .fna / .fasta / unknown → FASTA
    }
};

// Backward-compat alias: old code that used FastFastaReader still compiles.
using FastFastaReader = SeqReader;
