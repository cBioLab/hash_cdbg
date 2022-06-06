#ifndef COL_BOSS_F_ARRAY_HPP
#define COL_BOSS_F_ARRAY_HPP

#include <sdsl/int_vector.hpp>

using namespace sdsl;

class f_array_t {

private:
    int_vector<64> m_C = {0, 0, 0, 0, 0, 0, 0};

public:
    typedef int_vector<64>::size_type   size_type;
    typedef uint8_t                     symbol_type;
    const int_vector<64>& C =           m_C;

public:
    f_array_t() = default;
    explicit f_array_t(int_vector<64>& char_freqs):m_C(char_freqs){}
    size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;
    void load(std::istream& in);
    void swap(f_array_t& in);
    symbol_type operator[](size_type pos) const;
    size_type symbol_freq(symbol_type symbol) const;
    std::pair<symbol_type, size_type> inv_select(size_type pos) const;
    size_type symbol_select(symbol_type symbol, size_type rank) const;
};

inline f_array_t::size_type f_array_t::symbol_freq(f_array_t::symbol_type symbol) const {
    return C[symbol+1]-C[symbol];
}

inline f_array_t::symbol_type f_array_t::operator[](f_array_t::size_type pos) const {
    //TODO change this, it is not efficient
    if(pos<m_C[2]){
        return 1;
    }else if(m_C[2]<=pos && pos <m_C[3]){
        return 2;
    }else if(m_C[3]<=pos && pos <m_C[4]){
        return 3;
    }else if(m_C[4]<=pos && pos <m_C[5]){
        return 4;
    }else{
        return 5;
    }
}

inline std::pair<f_array_t::symbol_type, f_array_t::size_type> f_array_t::inv_select(f_array_t::size_type pos) const {
    symbol_type symbol = (*this)[pos];
    return {symbol, pos-C[symbol]+1};
}

inline f_array_t::size_type f_array_t::symbol_select(f_array_t::symbol_type symbol, f_array_t::size_type rank) const {
    return C[symbol] + rank;
}

#endif //COL_BOSS_F_ARRAY_HPP
