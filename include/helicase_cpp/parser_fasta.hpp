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

enum class Event { Record, DnaChunk };

// ─── FastaParser ─────────────────────────────────────────────────────────────
template <uint64_t CONFIG, typename I>
class FastaParser {
    FastaLexer<CONFIG, I> lexer_;

    // State machine
    enum class State { Start, Restart, Header, StartDNA, InDNABlock, EndDNA };
    State state_ = State::Start;

    bool finished_ = false;
    FastaChunk block_{};
    size_t block_counter_ = 0;
    size_t pos_in_block_ = 0;

    // Header
    size_t header_start_ = 0, header_end_ = 0;  // for RANDOM_ACCESS
    std::vector<uint8_t> cur_header_;

    // DNA
    size_t dna_start_ = 0, dna_end_ = 0;  // for RANDOM_ACCESS
    bool contiguous_dna_ = true;
    std::vector<uint8_t> cur_dna_string_;
    ColumnarDNA cur_dna_columnar_;
    PackedDNA cur_dna_packed_;
    BitMask cur_mask_non_actg_;
    BitMask cur_mask_n_;
    size_t dna_len_ = 0;

    size_t global_pos() const noexcept { return 64 * block_counter_ + pos_in_block_; }

    // ── Block advance helpers ─────────────────────────────────────────────
    bool advance_block() {
        FastaChunk c = lexer_.next_chunk();
        if (c.len == 0) return true;  // EOF
        block_ = c;
        block_counter_++;
        pos_in_block_ = 0;
        return false;
    }

    // Skip forward until `mask` has a set bit at or after pos_in_block_.
    // Returns true if EOF reached before finding a set bit.
    bool skip_to_mask(uint64_t FastaChunk::*mask_ptr) {
        uint64_t mask = (block_.*mask_ptr) & (~0ull << pos_in_block_);
        while (mask == 0) {
            if (advance_block()) return true;
            mask = block_.*mask_ptr;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        return false;
    }

    // ── skip_to_start_header ─────────────────────────────────────────────
    bool skip_to_start_header() {
        uint64_t mask = block_.header & (~0ull << pos_in_block_);
        while (mask == 0) {
            if (advance_block()) return true;
            mask = block_.header;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        return false;
    }

    // ── skip_to_header_or_dna ────────────────────────────────────────────
    bool skip_to_header_or_dna() {
        uint64_t mask = (block_.is_dna | block_.header) & (~0ull << pos_in_block_);
        while (mask == 0) {
            if (advance_block()) return true;
            mask = block_.is_dna | block_.header;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        return false;
    }

    // ── skip_to_end_header ───────────────────────────────────────────────
    bool skip_to_end_header() {
        uint64_t mask = ~block_.header & (~0ull << pos_in_block_);
        size_t first_pos = pos_in_block_;
        while (mask == 0) {
            if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && !I::RANDOM_ACCESS) {
                const uint8_t* blk = lexer_.input.current_block();
                cur_header_.insert(cur_header_.end(), blk + pos_in_block_, blk + block_.len);
            }
            if (advance_block()) return true;
            first_pos = 0;
            mask = ~block_.header;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && !I::RANDOM_ACCESS) {
            const uint8_t* blk = lexer_.input.current_block();
            cur_header_.insert(cur_header_.end(), blk + first_pos, blk + pos_in_block_);
        }
        return false;
    }

    // ── skip_to_non_dna ──────────────────────────────────────────────────
    // Accumulates DNA data while walking through is_dna bits.
    bool skip_to_non_dna() {
        uint64_t mask = ~block_.is_dna & (~0ull << pos_in_block_);
        size_t first_pos = pos_in_block_;
        while (mask == 0) {
            size_t chunk_size = block_.len - pos_in_block_;
            if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
                if (!I::RANDOM_ACCESS || !contiguous_dna_) {
                    const uint8_t* blk = lexer_.input.current_block();
                    cur_dna_string_.insert(cur_dna_string_.end(), blk + pos_in_block_, blk + block_.len);
                }
            }
            if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR)) {
                cur_dna_columnar_.append(block_.high_bit >> pos_in_block_, block_.low_bit >> pos_in_block_, chunk_size);
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
            if (advance_block()) {
                pos_in_block_ = block_.len;
                return true;
            }
            first_pos = 0;
            mask = ~block_.is_dna;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        size_t chunk_size = pos_in_block_ - first_pos;

        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
            if (!I::RANDOM_ACCESS || !contiguous_dna_) {
                const uint8_t* blk = lexer_.input.current_block();
                cur_dna_string_.insert(cur_dna_string_.end(), blk + first_pos, blk + pos_in_block_);
            }
        }
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR)) {
            cur_dna_columnar_.append(block_.high_bit >> first_pos, block_.low_bit >> first_pos, chunk_size);
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
        return false;
    }

    // ── skip_to_dna_or_split_or_header ──────────────────────────────────
    bool skip_to_dna_or_split_or_header() {
        uint64_t mask = (block_.is_dna | block_.split | block_.header) & (~0ull << pos_in_block_);
        while (mask == 0) {
            if (advance_block()) {
                pos_in_block_ = block_.len;
                return true;
            }
            mask = block_.is_dna | block_.split | block_.header;
        }
        pos_in_block_ = __builtin_ctzll(mask);
        return false;
    }

    // ── clear helpers ────────────────────────────────────────────────────
    void clear_record() {
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER))
            cur_header_.clear();
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
    // Construct from an input source.
    explicit FastaParser(I inp) : lexer_(std::move(inp)) {
        FastaChunk first = lexer_.next_chunk();
        if (first.len == 0) { finished_ = true; return; }
        block_ = first;
        block_counter_ = 0;
    }

    // Convenience: construct directly from (data, size) when I = SliceInput.
    FastaParser(const uint8_t* data, size_t size) : FastaParser(I(data, size)) {}
    FastaParser(const char* s, size_t size) : FastaParser(I((const uint8_t*)s, size)) {}
    explicit FastaParser(const char* s) : FastaParser(I((const uint8_t*)s, std::strlen(s))) {}

    // Convenience: construct from a byte slice.
    static FastaParser from_slice(const uint8_t* data, size_t size) {
        return FastaParser(I(data, size));
    }
    static FastaParser from_slice(const char* s) {
        return FastaParser(I((const uint8_t*)s, std::strlen(s)));
    }

    // ── Iterator ─────────────────────────────────────────────────────────
    // Returns the next event, or std::nullopt at end.
    std::optional<Event> next() {
        for (;;) {
            switch (state_) {
            case State::Start: {
                if (finished_) return std::nullopt;
                finished_ = skip_to_start_header();
                if (block_.header != 0) state_ = State::Header;
                break;
            }
            case State::Restart: {
                if (finished_) {
                    state_ = State::Start;
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::RETURN_RECORD))
                        return Event::Record;
                    continue;
                }
                finished_ = skip_to_header_or_dna();
                uint64_t bit = 1ull << pos_in_block_;
                if (bit & block_.header) {
                    state_ = State::Header;
                    if constexpr (advanced::flag_is_set(CONFIG, advanced::RETURN_RECORD))
                        return Event::Record;
                } else if (bit & block_.is_dna) {
                    state_ = State::StartDNA;
                }
                break;
            }
            case State::Header: {
                if constexpr (advanced::flag_is_not_set(CONFIG, advanced::MERGE_RECORDS))
                    clear_record();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && I::RANDOM_ACCESS)
                    header_start_ = global_pos() + 1;
                finished_ = skip_to_end_header();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_HEADER) && I::RANDOM_ACCESS)
                    header_end_ = global_pos() - 1;
                contiguous_dna_ = true;
                state_ = State::Restart;
                break;
            }
            case State::StartDNA: {
                state_ = State::InDNABlock;
                if constexpr (advanced::flag_is_not_set(CONFIG, advanced::MERGE_DNA_CHUNKS))
                    clear_chunk();
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING) && I::RANDOM_ACCESS)
                    dna_start_ = global_pos();
                break;
            }
            case State::InDNABlock: {
                if (skip_to_non_dna()) {
                    finished_ = true;
                    dna_end_ = global_pos();
                    state_ = State::EndDNA;
                    continue;
                }
                dna_end_ = global_pos();
                if (skip_to_dna_or_split_or_header()) {
                    finished_ = true;
                    state_ = State::EndDNA;
                    continue;
                }
                // If leaving DNA region and in RANDOM_ACCESS+contiguous mode, flush to string
                if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING) && I::RANDOM_ACCESS) {
                    if (contiguous_dna_ && !((1ull << pos_in_block_) & block_.header)) {
                        const uint8_t* data = lexer_.input.get_data();
                        cur_dna_string_.insert(cur_dna_string_.end(), data + dna_start_, data + dna_end_);
                        contiguous_dna_ = false;
                    }
                }
                if ((1ull << pos_in_block_) & block_.is_dna) {
                    state_ = State::InDNABlock;
                } else {
                    state_ = State::EndDNA;
                }
                break;
            }
            case State::EndDNA: {
                state_ = State::Restart;
                if constexpr (advanced::flag_is_set(CONFIG, advanced::RETURN_DNA_CHUNK))
                    return Event::DnaChunk;
                break;
            }
            }
        }
    }

    // ── Accessors ────────────────────────────────────────────────────────
    const std::vector<uint8_t>& get_header() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_HEADER))
            throw std::logic_error("headers are ignored");
        if constexpr (I::RANDOM_ACCESS) {
            // Return slice via range; for simplicity, use cur_header_ as staging buffer
            // Caller can also use get_header_range() for zero-copy.
        }
        return cur_header_;
    }

    // Zero-copy header access for RANDOM_ACCESS inputs.
    // Returns pointer + length into the underlying data.
    std::pair<const uint8_t*, size_t> get_header_raw() const {
        if constexpr (I::RANDOM_ACCESS) {
            return {lexer_.input.get_data() + header_start_, header_end_ - header_start_};
        }
        // For non-RA: skip leading '>' and trailing '\n'
        size_t n = cur_header_.size();
        size_t start = (n > 0 && cur_header_[0] == '>') ? 1 : 0;
        size_t end = (n > 0 && cur_header_[n-1] == '\n') ? n - 1 : n;
        return {cur_header_.data() + start, end - start};
    }

    std::vector<uint8_t> get_header_owned() {
        if constexpr (I::RANDOM_ACCESS) {
            auto [p, n] = get_header_raw();
            return std::vector<uint8_t>(p, p + n);
        }
        std::vector<uint8_t> res;
        std::swap(res, cur_header_);
        // strip leading '>' and trailing '\n'
        size_t start = (!res.empty() && res[0] == '>') ? 1 : 0;
        size_t end = (!res.empty() && res.back() == '\n') ? res.size() - 1 : res.size();
        return std::vector<uint8_t>(res.begin() + start, res.begin() + end);
    }

    const std::vector<uint8_t>& get_dna_string() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            throw std::logic_error("dna_string not enabled");
        return cur_dna_string_;
    }

    // Zero-copy DNA access: (ptr, len) valid until the next next() call.
    // For RANDOM_ACCESS + single-block DNA, points into the source buffer.
    // Otherwise points into the accumulated cur_dna_string_.
    std::pair<const char*, size_t> get_dna_raw() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            throw std::logic_error("dna_string not enabled");
        if constexpr (I::RANDOM_ACCESS) {
            if (contiguous_dna_ && dna_end_ > dna_start_)
                return {reinterpret_cast<const char*>(lexer_.input.get_data() + dna_start_),
                        dna_end_ - dna_start_};
        }
        return {reinterpret_cast<const char*>(cur_dna_string_.data()), cur_dna_string_.size()};
    }

    std::vector<uint8_t> get_dna_string_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_STRING))
            throw std::logic_error("dna_string not enabled");
        if constexpr (I::RANDOM_ACCESS) {
            if (contiguous_dna_) {
                const uint8_t* data = lexer_.input.get_data();
                return std::vector<uint8_t>(data + dna_start_, data + dna_end_);
            }
        }
        std::vector<uint8_t> res;
        std::swap(res, cur_dna_string_);
        return res;
    }

    const ColumnarDNA& get_dna_columnar() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR))
            throw std::logic_error("dna_columnar not enabled");
        return cur_dna_columnar_;
    }

    ColumnarDNA get_dna_columnar_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_COLUMNAR))
            throw std::logic_error("dna_columnar not enabled");
        ColumnarDNA res;
        std::swap(res, cur_dna_columnar_);
        return res;
    }

    const PackedDNA& get_dna_packed() const {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_PACKED))
            throw std::logic_error("dna_packed not enabled");
        return cur_dna_packed_;
    }

    PackedDNA get_dna_packed_owned() {
        if constexpr (advanced::flag_is_not_set(CONFIG, advanced::COMPUTE_DNA_PACKED))
            throw std::logic_error("dna_packed not enabled");
        PackedDNA res;
        std::swap(res, cur_dna_packed_);
        return res;
    }

    const BitMask& get_mask_non_actg() const { return cur_mask_non_actg_; }
    BitMask get_mask_non_actg_owned() {
        BitMask res; std::swap(res, cur_mask_non_actg_); return res;
    }

    const BitMask& get_mask_n() const { return cur_mask_n_; }
    BitMask get_mask_n_owned() {
        BitMask res; std::swap(res, cur_mask_n_); return res;
    }

    size_t get_dna_len() const {
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_LEN))
            return dna_len_;
        if constexpr (advanced::flag_is_set(CONFIG, advanced::COMPUTE_DNA_STRING)) {
            if constexpr (I::RANDOM_ACCESS) {
                if (contiguous_dna_) return dna_end_ - dna_start_;
            }
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
