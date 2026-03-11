#pragma once
#include <cstdint>

namespace helicase {

using Config = uint64_t;

namespace advanced {

[[nodiscard]] inline constexpr bool flag_is_set(Config cfg, Config flag) noexcept {
    return (cfg & flag) != 0;
}
[[nodiscard]] inline constexpr bool flag_is_not_set(Config cfg, Config flag) noexcept {
    return (cfg & flag) == 0;
}

static constexpr Config COMPUTE_HEADER        = 1ull << 0;
static constexpr Config COMPUTE_DNA_STRING    = 1ull << 1;
static constexpr Config COMPUTE_DNA_COLUMNAR  = 1ull << 2;
static constexpr Config COMPUTE_DNA_PACKED    = 1ull << 3;
static constexpr Config COMPUTE_DNA_LEN       = 1ull << 4;
static constexpr Config COMPUTE_QUALITY       = 1ull << 5;
static constexpr Config SPLIT_NON_ACTG        = 1ull << 6;
static constexpr Config RETURN_RECORD         = 1ull << 7;
static constexpr Config RETURN_DNA_CHUNK      = 1ull << 8;
static constexpr Config MERGE_DNA_CHUNKS      = 1ull << 9;
static constexpr Config MERGE_RECORDS         = 1ull << 10;
static constexpr Config COMPUTE_MASK_NON_ACTG = 1ull << 11;
static constexpr Config COMPUTE_MASK_N        = 1ull << 12;

static constexpr Config DEFAULT_CONFIG = COMPUTE_HEADER | COMPUTE_DNA_STRING | RETURN_RECORD;

} // namespace advanced

/// Compile-time builder for the parser configuration.
class ParserOptions {
    Config cfg_;

public:
    constexpr ParserOptions() noexcept : cfg_(advanced::DEFAULT_CONFIG) {}
    explicit constexpr ParserOptions(Config c) noexcept : cfg_(c) {}

    [[nodiscard]] constexpr Config config() const noexcept { return cfg_; }

    [[nodiscard]] constexpr ParserOptions compute_headers() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_HEADER);
    }
    [[nodiscard]] constexpr ParserOptions ignore_headers() const noexcept {
        return ParserOptions(cfg_ & ~advanced::COMPUTE_HEADER);
    }
    [[nodiscard]] constexpr ParserOptions compute_quality() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_QUALITY);
    }
    [[nodiscard]] constexpr ParserOptions ignore_quality() const noexcept {
        return ParserOptions(cfg_ & ~advanced::COMPUTE_QUALITY);
    }
    [[nodiscard]] constexpr ParserOptions compute_dna_len() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_DNA_LEN);
    }
    [[nodiscard]] constexpr ParserOptions ignore_dna() const noexcept {
        return ParserOptions(cfg_ & ~(advanced::COMPUTE_DNA_STRING | advanced::COMPUTE_DNA_COLUMNAR |
                                      advanced::COMPUTE_DNA_PACKED | advanced::SPLIT_NON_ACTG |
                                      advanced::RETURN_DNA_CHUNK));
    }
    [[nodiscard]] constexpr ParserOptions dna_string() const noexcept {
        return ignore_dna().and_dna_string();
    }
    [[nodiscard]] constexpr ParserOptions and_dna_string() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_DNA_STRING);
    }
    [[nodiscard]] constexpr ParserOptions dna_packed() const noexcept {
        return ignore_dna().and_dna_packed();
    }
    [[nodiscard]] constexpr ParserOptions and_dna_packed() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_DNA_PACKED | advanced::SPLIT_NON_ACTG |
                              advanced::RETURN_DNA_CHUNK);
    }
    [[nodiscard]] constexpr ParserOptions dna_columnar() const noexcept {
        return ignore_dna().and_dna_columnar();
    }
    [[nodiscard]] constexpr ParserOptions and_dna_columnar() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_DNA_COLUMNAR | advanced::SPLIT_NON_ACTG |
                              advanced::RETURN_DNA_CHUNK);
    }
    [[nodiscard]] constexpr ParserOptions compute_mask_non_actg() const noexcept {
        return ParserOptions(cfg_ | advanced::COMPUTE_MASK_NON_ACTG);
    }
    [[nodiscard]] constexpr ParserOptions keep_non_actg() const noexcept {
        return ParserOptions(cfg_ & ~(advanced::SPLIT_NON_ACTG | advanced::RETURN_DNA_CHUNK |
                                       advanced::MERGE_DNA_CHUNKS));
    }
    [[nodiscard]] constexpr ParserOptions split_non_actg() const noexcept {
        return ParserOptions((cfg_ & ~advanced::MERGE_DNA_CHUNKS) | advanced::SPLIT_NON_ACTG |
                              advanced::RETURN_DNA_CHUNK);
    }
    [[nodiscard]] constexpr ParserOptions skip_non_actg() const noexcept {
        return ParserOptions((cfg_ & ~advanced::RETURN_DNA_CHUNK) | advanced::SPLIT_NON_ACTG |
                              advanced::MERGE_DNA_CHUNKS);
    }
    [[nodiscard]] constexpr ParserOptions return_record(bool enable) const noexcept {
        return enable ? ParserOptions(cfg_ | advanced::RETURN_RECORD)
                      : ParserOptions(cfg_ & ~advanced::RETURN_RECORD);
    }
    [[nodiscard]] constexpr ParserOptions return_dna_chunk(bool enable) const noexcept {
        return enable ? ParserOptions(cfg_ | advanced::RETURN_DNA_CHUNK)
                      : ParserOptions(cfg_ & ~advanced::RETURN_DNA_CHUNK);
    }

    // operator& for combining with additional flag masks (e.g. & ~RETURN_RECORD)
    [[nodiscard]] constexpr ParserOptions operator&(Config mask) const noexcept {
        return ParserOptions(cfg_ & mask);
    }
    [[nodiscard]] constexpr ParserOptions operator|(Config mask) const noexcept {
        return ParserOptions(cfg_ | mask);
    }
    [[nodiscard]] constexpr Config operator~() const noexcept { return ~cfg_; }
};

} // namespace helicase
