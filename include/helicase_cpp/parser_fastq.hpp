#pragma once
#include "config.hpp"
#include "dna_format.hpp"
#include "input.hpp"
#include "lexer.hpp"
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace helicase {

// ─── FastqParser ─────────────────────────────────────────────────────────────
template <uint64_t CONFIG, typename I>
class FastqParser {
    FastqLexer<CONFIG, I> lexer_;

    bool finished_ = false;
    int line_count_ = 0;
    FastqChunk block_{};
    size_t block_counter_ = 0;
    size_t pos_in_block_ = 0;

    // Header
    size_t header_start_ = 0, header_end_ = 0;
    std::vector<uint8_t> cur_header_;

    // Quality
    size_t quality_start_ = 0, quality_end_ = 0;
    std::vector<uint8_t> cur_quality_;

    // DNA
    size_t dna_start_ = 0, dna_end_ = 0;
    std::vector<uint8_t> cur_dna_string_;
    ColumnarDNA cur_dna_columnar_;
    PackedDNA cur_dna_packed_;
    BitMask cur_mask_non_actg_;
    BitMask cur_mask_n_;
    size_t dna_len_ = 0;

    size_t global_pos() const noexcept { return 64 * block_counter_ + pos_in_block_; }

    // Advance to the next block.
    bool advance_block() {
        FastqChunk c = lexer_.next_chunk();
        if (c.len == 0) { finished_ = true; return true; }
        block_ = c;
        block_counter_++;
        pos_in_block_ = 0;
        return false;
    }

    // Move pos_in_block_ forward by one (or advance to next block).
    void increment_pos() {
        if (pos_in_block_ + 1 < block_.len) {
            pos_in_block_++;
        } else {
            advance_block();
        }
    }

    // Consume the lowest set newline bit, advance pos, increment line_count.
    void consume_newline() {
        block_.newline &= block_.newline - 1;  // clear lowest set bit
        increment_pos();
        line_count_++;
    }

    void clear_record() {
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER))
            cur_header_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY))
            cur_quality_.clear();
        clear_chunk();
    }

    void clear_chunk() {
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            cur_dna_string_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR))
            cur_dna_columnar_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_PACKED))
            cur_dna_packed_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_NON_ACTG))
            cur_mask_non_actg_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_N))
            cur_mask_n_.clear();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_LEN))
            dna_len_ = 0;
    }

public:
    explicit FastqParser(I inp) : lexer_(std::move(inp)) {
        FastqChunk first = lexer_.next_chunk();
        if (first.len == 0) { finished_ = true; return; }
        block_ = first;
        block_counter_ = 0;
    }

    FastqParser(const uint8_t* data, size_t size) : FastqParser(I(data, size)) {}
    FastqParser(const char* s, size_t size) : FastqParser(I((const uint8_t*)s, size)) {}
    explicit FastqParser(const char* s) : FastqParser(I((const uint8_t*)s, std::strlen(s))) {}

    static FastqParser from_slice(const uint8_t* data, size_t size) {
        return FastqParser(I(data, size));
    }
    static FastqParser from_slice(const char* s) {
        return FastqParser(I((const uint8_t*)s, std::strlen(s)));
    }

    std::optional<Event> next() {
        for (;;) {
            switch (line_count_ % 4) {
            case 0: {
                // HEADER LINE: skip '@', collect until newline
                increment_pos();
                if (finished_) return std::nullopt;
                if constexpr (advanced::flag_is_not_set(CONFIG, advanced::MERGE_RECORDS))
                    clear_record();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && I::RANDOM_ACCESS)
                    header_start_ = global_pos();
                size_t first_pos = pos_in_block_;
                while (block_.newline == 0) {
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && !I::RANDOM_ACCESS) {
                        const uint8_t* blk = lexer_.input.current_block();
                        cur_header_.insert(cur_header_.end(), blk + first_pos, blk + block_.len);
                    }
                    if (advance_block()) return std::nullopt;
                    first_pos = 0;
                }
                pos_in_block_ = __builtin_ctzll(block_.newline);
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && I::RANDOM_ACCESS) {
                    header_end_ = global_pos();
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && !I::RANDOM_ACCESS) {
                    const uint8_t* blk = lexer_.input.current_block();
                    cur_header_.insert(cur_header_.end(), blk + first_pos, blk + pos_in_block_);
                }
                consume_newline();
                break;
            }
            case 1: {
                // SEQUENCE LINE
                if constexpr (advanced::flag_is_set(CONFIG, advanced::SPLIT_NON_ACTG)) {
                    // Skip to first DNA or newline
                    uint64_t mask = (block_.is_dna | block_.newline) & (~0ull << pos_in_block_);
                    while (mask == 0) {
                        if (advance_block()) return std::nullopt;
                        mask = block_.is_dna | block_.newline;
                    }
                    pos_in_block_ = __builtin_ctzll(mask);
                    if ((1ull << pos_in_block_) & block_.newline) {
                        consume_newline();
                        continue;
                    }
                }

                // Collect DNA until non-DNA position
                uint64_t mask = ~block_.is_dna & (~0ull << pos_in_block_);
                if constexpr (advanced::flag_is_not_set(CONFIG, advanced::MERGE_DNA_CHUNKS))
                    clear_chunk();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING) &&
                              advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG) && I::RANDOM_ACCESS) {
                    dna_start_ = global_pos();
                }
                size_t first_pos = pos_in_block_;
                while (mask == 0) {
                    size_t chunk_size = block_.len - pos_in_block_;
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
                        if constexpr (advanced::flag_is_set(CONFIG, advanced::SPLIT_NON_ACTG) || !I::RANDOM_ACCESS) {
                            const uint8_t* blk = lexer_.input.current_block();
                            cur_dna_string_.insert(cur_dna_string_.end(), blk + pos_in_block_, blk + block_.len);
                        }
                    }
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR)) {
                        cur_dna_columnar_.append(block_.high_bit >> pos_in_block_,
                                                  block_.low_bit  >> pos_in_block_, chunk_size);
                    }
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_PACKED)) {
                        cur_dna_packed_.append(block_.two_bits >> (2 * pos_in_block_), 2 * chunk_size);
                    }
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_NON_ACTG)) {
                        cur_mask_non_actg_.append(block_.mask_non_actg >> pos_in_block_, chunk_size);
                    }
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_N)) {
                        cur_mask_n_.append(block_.mask_n >> pos_in_block_, chunk_size);
                    }
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_LEN)) {
                        dna_len_ += chunk_size;
                    }
                    if (advance_block()) return std::nullopt;
                    first_pos = 0;
                    mask = ~block_.is_dna;
                }
                pos_in_block_ = __builtin_ctzll(mask);
                size_t chunk_size = pos_in_block_ - first_pos;

                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::SPLIT_NON_ACTG) || !I::RANDOM_ACCESS) {
                        const uint8_t* blk = lexer_.input.current_block();
                        cur_dna_string_.insert(cur_dna_string_.end(), blk + first_pos, blk + pos_in_block_);
                    }
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR)) {
                    cur_dna_columnar_.append(block_.high_bit >> first_pos,
                                              block_.low_bit  >> first_pos, chunk_size);
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_PACKED)) {
                    cur_dna_packed_.append(block_.two_bits >> (2 * first_pos), 2 * chunk_size);
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_NON_ACTG)) {
                    cur_mask_non_actg_.append(block_.mask_non_actg >> first_pos, chunk_size);
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_MASK_N)) {
                    cur_mask_n_.append(block_.mask_n >> first_pos, chunk_size);
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_LEN)) {
                    dna_len_ += chunk_size;
                }

                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING) &&
                              advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG) && I::RANDOM_ACCESS) {
                    dna_end_ = global_pos();
                }
                // Consume newline if we're at it
                if (advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG) ||
                    ((1ull << pos_in_block_) & block_.newline)) {
                    consume_newline();
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::RETURN_DNA_CHUNK))
                    return Event::DnaChunk;
                break;
            }
            case 2: {
                // PLUS LINE: skip to newline
                while (block_.newline == 0) {
                    if (advance_block()) return std::nullopt;
                }
                pos_in_block_ = __builtin_ctzll(block_.newline);
                consume_newline();
                break;
            }
            case 3: {
                // QUALITY LINE
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY) && I::RANDOM_ACCESS)
                    quality_start_ = global_pos();
                size_t first_pos = pos_in_block_;
                bool eof_in_quality = false;
                while (block_.newline == 0) {
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY) && !I::RANDOM_ACCESS) {
                        const uint8_t* blk = lexer_.input.current_block();
                        cur_quality_.insert(cur_quality_.end(), blk + pos_in_block_, blk + block_.len);
                    }
                    if (advance_block()) {
                        eof_in_quality = true;
                        break;  // EOF: fall through to finalize
                    }
                    first_pos = 0;
                }
                // After loop: block_.newline may be 0 (EOF) or nonzero (found newline).
                // Match Rust behavior: trailing_zeros() of 0 = 64 on u64, clamped to block_.len.
                if (block_.newline != 0) {
                    pos_in_block_ = __builtin_ctzll(block_.newline);
                } else {
                    pos_in_block_ = 64;  // will be clamped below
                }
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY))
                    pos_in_block_ = std::min(pos_in_block_, block_.len);
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY) && I::RANDOM_ACCESS)
                    quality_end_ = global_pos();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_QUALITY) && !I::RANDOM_ACCESS) {
                    if (!eof_in_quality) {
                        const uint8_t* blk = lexer_.input.current_block();
                        cur_quality_.insert(cur_quality_.end(), blk + first_pos, blk + pos_in_block_);
                    }
                }
                consume_newline();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::RETURN_RECORD))
                    return Event::Record;
                break;
            }
            }
        }
    }

    // ── Accessors ────────────────────────────────────────────────────────
    std::vector<uint8_t> get_header_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_HEADER))
            throw std::logic_error("headers are ignored");
        if constexpr (I::RANDOM_ACCESS) {
            const uint8_t* data = lexer_.input.get_data();
            // Skip the leading '@' character
            size_t start = header_start_;
            if (start < lexer_.input.get_size() && data[start] == '@') start++;
            return std::vector<uint8_t>(data + start, data + header_end_);
        }
        std::vector<uint8_t> res;
        std::swap(res, cur_header_);
        // Strip leading '@' if present
        size_t start = (!res.empty() && res[0] == '@') ? 1 : 0;
        return std::vector<uint8_t>(res.begin() + start, res.end());
    }

    const std::vector<uint8_t>& get_quality_ref() const { return cur_quality_; }

    std::optional<std::vector<uint8_t>> get_quality_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_QUALITY))
            throw std::logic_error("quality is ignored");
        if constexpr (I::RANDOM_ACCESS) {
            const uint8_t* data = lexer_.input.get_data();
            return std::vector<uint8_t>(data + quality_start_, data + quality_end_);
        }
        std::vector<uint8_t> res;
        std::swap(res, cur_quality_);
        return res;
    }

    // Zero-copy DNA access: (ptr, len) valid until the next next() call.
    // For RANDOM_ACCESS + non-SPLIT_NON_ACTG (single-line FASTQ), points into the source buffer.
    // Otherwise points into the accumulated cur_dna_string_.
    std::pair<const char*, size_t> get_dna_raw() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            throw std::logic_error("dna_string not enabled");
        if constexpr (I::RANDOM_ACCESS && advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG)) {
            if (dna_end_ > dna_start_)
                return {reinterpret_cast<const char*>(lexer_.input.get_data() + dna_start_),
                        dna_end_ - dna_start_};
        }
        return {reinterpret_cast<const char*>(cur_dna_string_.data()), cur_dna_string_.size()};
    }

    std::vector<uint8_t> get_dna_string_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            throw std::logic_error("dna_string not enabled");
        if constexpr (I::RANDOM_ACCESS && advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG)) {
            const uint8_t* data = lexer_.input.get_data();
            return std::vector<uint8_t>(data + dna_start_, data + dna_end_);
        }
        std::vector<uint8_t> res;
        std::swap(res, cur_dna_string_);
        return res;
    }

    ColumnarDNA get_dna_columnar_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR))
            throw std::logic_error("dna_columnar not enabled");
        ColumnarDNA res; std::swap(res, cur_dna_columnar_); return res;
    }

    PackedDNA get_dna_packed_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_PACKED))
            throw std::logic_error("dna_packed not enabled");
        PackedDNA res; std::swap(res, cur_dna_packed_); return res;
    }

    size_t get_dna_len() const {
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_LEN)) return dna_len_;
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
            if constexpr (I::RANDOM_ACCESS && advanced::flag_is_not_set(CONFIG, advanced::SPLIT_NON_ACTG))
                return dna_end_ - dna_start_;
            return cur_dna_string_.size();
        }
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR))
            return cur_dna_columnar_.len();
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_PACKED))
            return cur_dna_packed_.len();
        throw std::logic_error("dna is ignored");
    }
};

} // namespace helicase
