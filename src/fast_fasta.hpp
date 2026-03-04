#pragma once

// SeqReader — unified FASTA/FASTQ reader with optional gzip decompression.
//
// Drop-in replacement for the old FastFastaReader (now supports FASTQ and .gz).
//
// Design: one instance per worker thread, kept alive across files.  call load()
// for each new file; the internal buffers only grow, eliminating per-file
// allocation overhead.  Sequences are returned via seq() / seq_len() and are
// valid until the next call to read_next_seq() or load().
//
// Format detection (from file extension, case-insensitive; .gz stripped first):
//   FASTA : .fa  .fna  .fasta  (and any other extension)
//   FASTQ : .fq  .fastq
//
// Gzip decompression: zlib gzFile — the compressed file is fully decompressed
// into buf_ on load(), then parsed in place.  This keeps the hot parse loop
// identical to the uncompressed path.
//
// FASTQ: single-line sequence format (standard illumina output).  Quality lines
// are skipped with a single memchr each — no per-base overhead.

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

    // Load a new file.  Auto-detects format and compression from extension.
    // Reuses existing buffer capacity.  Returns false if the file cannot be
    // opened or yields no data.
    bool load(const std::string& path)
    {
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
        seq_.clear();
        return fmt_ == SeqFormat::FASTQ ? read_fastq() : read_fasta();
    }

    const char* seq()     const { return seq_.data(); }
    size_t      seq_len() const { return seq_.size(); }
    void        close()         {} // no-op; kept for drop-in compatibility

private:

    SeqFormat         fmt_ = SeqFormat::FASTA;
    std::vector<char> buf_;   // raw (or decompressed) file content; reused
    std::vector<char> seq_;   // current stripped sequence; reused
    size_t            rpos_ = 0; // read cursor in buf_


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
        gzFile f = gzopen(path.c_str(), "rb");
        if (!f) return false;
        gzbuffer(f, 256u << 10); // 256 KB zlib internal buffer

        // Decompress into buf_ using a doubling-resize strategy.
        buf_.resize(1u << 20); // start at 1 MB
        size_t total = 0;
        int n;
        while ((n = gzread(f, buf_.data() + total,
                           static_cast<unsigned>(buf_.size() - total))) > 0) {
            total += static_cast<size_t>(n);
            if (total == buf_.size())
                buf_.resize(buf_.size() * 2);
        }
        gzclose(f);
        buf_.resize(total);
        if (total == 0) return false;
        rpos_ = 0;
        return true;
    }


    // ── Parsers ───────────────────────────────────────────────────────────────

    // FASTA: collect sequence lines between '>' headers, stripping newlines.
    bool read_fasta()
    {
        // Seek to the next '>'
        while (rpos_ < buf_.size() && buf_[rpos_] != '>') ++rpos_;
        if (rpos_ >= buf_.size()) return false;
        ++rpos_; // consume '>'

        // Skip header line
        const char* nl = static_cast<const char*>(
            std::memchr(buf_.data() + rpos_, '\n', buf_.size() - rpos_));
        if (!nl) return false;
        rpos_ = static_cast<size_t>(nl - buf_.data()) + 1;

        // Collect sequence bytes line by line until the next '>' or EOF.
        while (rpos_ < buf_.size() && buf_[rpos_] != '>') {
            const char*  start     = buf_.data() + rpos_;
            const size_t remaining = buf_.size() - rpos_;

            nl = static_cast<const char*>(std::memchr(start, '\n', remaining));
            size_t line_len = nl ? static_cast<size_t>(nl - start) : remaining;
            if (line_len > 0 && start[line_len - 1] == '\r') --line_len;

            if (line_len > 0) {
                const size_t old_sz = seq_.size();
                seq_.resize(old_sz + line_len);
                std::memcpy(seq_.data() + old_sz, start, line_len);
            }

            rpos_ += nl ? static_cast<size_t>(nl - start) + 1 : remaining;
        }

        return !seq_.empty();
    }

    // FASTQ: 4-line records — @header / sequence / +ignored / quality.
    // Only the sequence line is kept; quality is skipped with a single memchr.
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
