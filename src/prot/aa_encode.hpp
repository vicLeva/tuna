#pragma once
#include <array>
#include <cstdint>

namespace prot {

// 20 standard amino acids in alphabetical order: A C D E F G H I K L M N P Q R S T V W Y
// Encoded 0-19. All others (B, J, O, U, X, Z, *, -) → AA_INVALID (sequence delimiter).

static constexpr uint8_t AA_INVALID = 0xFF;
static constexpr uint8_t AA_COUNT   = 20;

inline constexpr std::array<uint8_t, 128> make_encode_table() noexcept {
    std::array<uint8_t, 128> t{};
    for (auto& v : t) v = AA_INVALID;
    // Standard 20 amino acids, uppercase
    t['A']=0;  t['C']=1;  t['D']=2;  t['E']=3;  t['F']=4;
    t['G']=5;  t['H']=6;  t['I']=7;  t['K']=8;  t['L']=9;
    t['M']=10; t['N']=11; t['P']=12; t['Q']=13; t['R']=14;
    t['S']=15; t['T']=16; t['V']=17; t['W']=18; t['Y']=19;
    // Lowercase (same encoding)
    t['a']=0;  t['c']=1;  t['d']=2;  t['e']=3;  t['f']=4;
    t['g']=5;  t['h']=6;  t['i']=7;  t['k']=8;  t['l']=9;
    t['m']=10; t['n']=11; t['p']=12; t['q']=13; t['r']=14;
    t['s']=15; t['t']=16; t['v']=17; t['w']=18; t['y']=19;
    return t;
}

static constexpr auto AA_ENCODE = make_encode_table();

static constexpr char AA_DECODE[20] = {
    'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'
};

inline bool is_aa(char c) noexcept {
    return static_cast<uint8_t>(c) < 128 && AA_ENCODE[static_cast<uint8_t>(c)] != AA_INVALID;
}

inline uint8_t encode_aa(char c) noexcept {
    return AA_ENCODE[static_cast<uint8_t>(c) & 0x7Fu];
}

// Decode a protein k-mer (5-bit packed uint64_t) to a string buffer.
// buf must have space for k+1 chars. Writes a NUL terminator.
template <uint16_t k>
inline void decode_pkmer(uint64_t kmer, char* buf) noexcept {
    for (int i = k - 1; i >= 0; --i) {
        buf[i] = AA_DECODE[kmer & 0x1fu];
        kmer >>= 5;
    }
    buf[k] = '\0';
}

} // namespace prot
