#ifndef COL_BOSS_DNA_ALPHABET_HPP
#define COL_BOSS_DNA_ALPHABET_HPP

#include <sdsl/int_vector.hpp>

class dna_alphabet {
public:

    typedef sdsl::int_vector<8>        char_arr_type;
    typedef uint8_t                    symbol_type;

    static const size_t                sigma=6; //{#,$,A,C,G,arr_t} //TODO including the N symbol is pending
    static const char_arr_type         comp2char;
    static const char_arr_type         comp2rev;
    static const char_arr_type         char2comp;
    static const char_arr_type         char2rev;

    dna_alphabet() = default;
};


#endif //COL_BOSS_DNA_ALPHABET_HPP
