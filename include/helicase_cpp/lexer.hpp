#pragma once
#include "carrying_add.hpp"
#include "config.hpp"
#include "input.hpp"
#include "simd.hpp"

namespace helicase {

// ─── FastaChunk ──────────────────────────────────────────────────────────────
struct FastaChunk {
    size_t len = 0;
    uint64_t header = 0;
    uint64_t split = 0;
    uint64_t is_dna = 0;
    __uint128_t two_bits = 0;
    uint64_t high_bit = 0;
    uint64_t low_bit = 0;
    uint64_t mask_non_actg = 0;
    uint64_t mask_n = 0;

    bool operator==(const FastaChunk&) const = default;
};

// ─── FastqChunk ──────────────────────────────────────────────────────────────
struct FastqChunk {
    size_t len = 0;
    uint64_t newline = 0;
    uint64_t is_dna = 0;
    __uint128_t two_bits = 0;
    uint64_t high_bit = 0;
    uint64_t low_bit = 0;
    uint64_t mask_non_actg = 0;
    uint64_t mask_n = 0;
};

// ─── FastaLexer ──────────────────────────────────────────────────────────────
// Produces FastaChunks by extracting bitmasks from 64-byte blocks.
template <uint64_t CONFIG, typename I>
class FastaLexer {
public:
    I input;

    explicit FastaLexer(I inp) : input(std::move(inp)) {}

    // Returns the next chunk, or a default (zero) FastaChunk with len==0 at EOF.
    FastaChunk next_chunk() {
        Block blk = input.next();
        if (!blk) return {};
        // Zero-pad to 64 bytes for SIMD (already done by input for last block)
        alignas(64) uint8_t buf[64] = {};
        std::memcpy(buf, blk.ptr, blk.size);

        FastaBitmask m = simd::extract_fasta_bitmask<CONFIG>(buf);
        uint64_t non_lf = ~m.line_feeds;
        uint64_t c = carry_.add(m.open_bracket, non_lf);
        uint64_t header = c ^ non_lf;
        uint64_t is_dna = m.is_dna & ~header & non_lf;
        uint64_t split = 0;
        if constexpr (advanced::flag_is_set(CONFIG, advanced::SPLIT_NON_ACTG))
            split = ~header & ~is_dna & non_lf;

        return FastaChunk{
            blk.size, header, split, is_dna,
            m.two_bits, m.high_bit, m.low_bit,
            m.mask_non_actg & ~m.line_feeds,
            m.mask_n
        };
    }

private:
    Carry carry_{false};
};

// ─── FastqLexer ──────────────────────────────────────────────────────────────
// Produces FastqChunks by extracting bitmasks from 64-byte blocks.
template <uint64_t CONFIG, typename I>
class FastqLexer {
public:
    I input;

    explicit FastqLexer(I inp) : input(std::move(inp)) {}

    FastqChunk next_chunk() {
        Block blk = input.next();
        if (!blk) return {};
        alignas(64) uint8_t buf[64] = {};
        std::memcpy(buf, blk.ptr, blk.size);

        FastqBitmask m = simd::extract_fastq_bitmask<CONFIG>(buf);
        uint64_t is_dna = m.is_dna & ~m.line_feeds;

        return FastqChunk{
            blk.size,
            m.line_feeds,
            is_dna,
            m.two_bits, m.high_bit, m.low_bit,
            m.mask_non_actg & ~m.line_feeds,
            m.mask_n
        };
    }
};

} // namespace helicase
