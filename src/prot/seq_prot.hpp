#pragma once

// ProtSeqSource — streaming protein FASTA reader (plain or gzip).
//
// Reads FASTA files, skips '>' header lines, calls on_chunk(const char*, size_t)
// for every run of valid amino acid characters within each sequence.
// Splits on: sequence boundaries ('>'), ambiguous residues (B,J,O,U,X,Z,*,-).
// Newlines inside sequence lines are ignored (multi-line FASTA support).

#include "aa_encode.hpp"
#include <zlib.h>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace prot {

// ── split_protein ─────────────────────────────────────────────────────────────

template <typename F>
inline void split_protein(const char* seq, size_t len, F&& on_chunk) {
    size_t start  = 0;
    bool   in_run = false;
    for (size_t i = 0; i < len; ++i) {
        if (is_aa(seq[i])) {
            if (!in_run) { start = i; in_run = true; }
        } else {
            if (in_run) { on_chunk(seq + start, i - start); in_run = false; }
        }
    }
    if (in_run) on_chunk(seq + start, len - start);
}


// ── ProtSeqSource ─────────────────────────────────────────────────────────────
//
// State machine:
//   in_header_  — consuming a header line (skip until '\n')
//   carry_      — partial AA run carried across block boundaries
//                 (flushed and emitted at: '>', non-AA delimiter, end-of-sequence)
//
// Invariant: carry_ always holds a contiguous run of valid AAs from the
// current sequence.  It is flushed (emitted) before any sequence boundary.

struct ProtSeqSource {
    static constexpr size_t BUFSIZE = 1u << 18;  // 256 KB

    ProtSeqSource() = default;

    template <typename F>
    void process(const std::string& path, F&& on_chunk) {
        // Reset per-file state.
        in_header_ = false;
        carry_.clear();

        const bool is_gz = path.size() > 3 &&
                           path.compare(path.size() - 3, 3, ".gz") == 0;
        if (is_gz)
            process_gz(path, std::forward<F>(on_chunk));
        else
            process_plain(path, std::forward<F>(on_chunk));

        // Flush any trailing carry at end of file.
        flush_carry(on_chunk);
        carry_.clear();
    }

private:
    bool        in_header_ = false;
    std::string carry_;

    // Flush and emit carry_ if non-empty.
    template <typename F>
    void flush_carry(F&& on_chunk) {
        if (!carry_.empty()) {
            on_chunk(carry_.data(), carry_.size());
            carry_.clear();
        }
    }

    // Process one raw buffer block.
    template <typename F>
    void parse_block(const char* buf, size_t len, F&& on_chunk) {
        size_t run_start = 0;
        bool   in_run    = false;

        auto flush_run = [&](size_t end) {
            if (!in_run) return;
            const char* p   = buf + run_start;
            size_t      rlen = end - run_start;
            if (!carry_.empty()) {
                carry_.append(p, rlen);
                on_chunk(carry_.data(), carry_.size());
                carry_.clear();
            } else {
                on_chunk(p, rlen);
            }
            in_run = false;
        };

        auto break_sequence = [&](size_t i) {
            // End current in-block run.
            flush_run(i);
            // Discard any cross-block carry (sequence boundary or invalid AA).
            carry_.clear();
        };

        for (size_t i = 0; i < len; ++i) {
            const char c = buf[i];

            if (in_header_) {
                break_sequence(i);   // ensure nothing leaks from before header
                if (c == '\n') in_header_ = false;
                continue;
            }

            if (c == '>') {
                break_sequence(i);
                in_header_ = true;
                continue;
            }

            if (c == '\n' || c == '\r') {
                // Newlines don't break the sequence run; just accumulate to carry.
                if (in_run) {
                    carry_.append(buf + run_start, i - run_start);
                    in_run = false;
                }
                continue;
            }

            if (is_aa(c)) {
                if (!in_run) { run_start = i; in_run = true; }
            } else {
                // Non-AA delimiter (B, X, *, -, etc.).
                flush_run(i);
                carry_.clear();
            }
        }

        // End of block: save in-block run tail to carry_ for next block.
        if (in_run) {
            carry_.append(buf + run_start, len - run_start);
        }
    }

    template <typename F>
    void process_plain(const std::string& path, F&& on_chunk) {
        FILE* fp = fopen(path.c_str(), "rb");
        if (!fp) throw std::runtime_error("Cannot open: " + path);
        buf_.resize(BUFSIZE);
        size_t n;
        while ((n = fread(buf_.data(), 1, BUFSIZE, fp)) > 0)
            parse_block(buf_.data(), n, on_chunk);
        fclose(fp);
    }

    template <typename F>
    void process_gz(const std::string& path, F&& on_chunk) {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (!gz) throw std::runtime_error("Cannot open gz: " + path);
        gzbuffer(gz, 256u << 10);
        buf_.resize(BUFSIZE);
        int n;
        while ((n = gzread(gz, buf_.data(), static_cast<unsigned>(BUFSIZE))) > 0)
            parse_block(buf_.data(), static_cast<size_t>(n), on_chunk);
        gzclose(gz);
    }

    std::vector<char> buf_;
};

} // namespace prot
