#ifndef COL_BOSS_HPP
#define COL_BOSS_HPP

#include <iostream>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/wt_huff.hpp>
#include <sdsl/vlc_vector.hpp>
#include "f_array_t.hpp"
#include "dna_alphabet.hpp"
#include "color_matrix.hpp"

class dbg_boss{

public:
    typedef uint64_t                                                   size_type;
    typedef std::pair<size_type, size_type>                            range_t;
    typedef std::vector<uint8_t>                                       label_t;
    typedef color_matrix<sdsl::int_vector<64>, sdsl::rrr_vector<63>, dbg_boss>   color_matrix_t;  // 16 to 64

    //degree result
    struct degree_t{
        size_type degree; //number of incomming/outgoing edges
        bool dollar_edge; //true if the dBG has incomming or outgoing edge labelled with $
    };

    //backward search result
    struct bs_t{
        size_type r_start; //start of the range
        size_type r_end; //end of the range
        size_type mm; //number of mismatches
        bs_t(size_type start, size_type end, size_type mis): r_start(start),
                                                             r_end(end),
                                                             mm(mis){}
    };


private:
    typedef rrr_vector<63>                              comp_bv_t;
    typedef typename sdsl::wt_huff<>                    edges_t;
    typedef typename rrr_vector<63>::rank_0_type        rank_0_t;
    typedef typename rrr_vector<63>::rank_1_type        rank_1_t;
    typedef typename rrr_vector<63>::select_0_type      select_0_t;
    typedef typename rrr_vector<63>::select_1_type      select_1_t;

    //backward search step
    struct bs_step_t{
        size_type r_start;
        size_type r_end;
        size_type mm;
        size_type last_index;
    };

    size_type                                           m_n_solid_nodes;
    size_t                                              m_k;
    f_array_t                                           m_f_array;
    edges_t                                             m_edge_bwt;
    comp_bv_t                                           m_node_marks;
    comp_bv_t                                           m_solid_nodes;
    rank_0_t                                            m_node_marks_rs;
    select_0_t                                          m_node_marks_ss;
    select_1_t                                          m_solid_nodes_ss;
    rank_1_t                                            m_solid_nodes_rs;
    const uint8_t                                       bit_clear=254;

    //get_colors data
    color_matrix_t                                      m_color_matrix;

    //temporal variables for interval_symbols function
    static std::vector<edges_t::value_type> cs;
    static std::vector<edges_t::size_type> c_i;
    static std::vector<edges_t::size_type> c_j;

public:
    const size_t&       k                 =             m_k;
    const size_type&    n_solid_nodes     =             m_n_solid_nodes;
    const f_array_t&    f_array           =             m_f_array;
    const edges_t&      edge_bwt          =             m_edge_bwt;
    const comp_bv_t&    node_marks        =             m_node_marks;
    const rank_0_t&     node_marks_rs     =             m_node_marks_rs;
    const select_0_t&   node_marks_ss     =             m_node_marks_ss;
    const comp_bv_t&    solid_nodes       =             m_solid_nodes;
    const select_1_t&   solid_nodes_ss    =             m_solid_nodes_ss;
    const rank_1_t&     solid_nodes_rs    =             m_solid_nodes_rs;

public:
    //constructors
    dbg_boss();
    dbg_boss(std::string &input_file, size_t K, size_t n_threads);
    dbg_boss(size_t K, size_t n_threads, std::string dir, std::string name);

    //navigational functions
    size_type outgoing(size_type v, uint8_t value, bool is_symbol=false,
                       std::vector<edges_t::value_type> &syms=cs,
                       std::vector<edges_t::size_type> &r_b=c_i,
                       std::vector<edges_t::size_type> &r_a=c_j) const;

    size_type incomming(size_type v, size_t value, bool is_symbol=false) const;
    size_t indegree(size_type v) const;
    size_t outdegree(size_type v) const;
    std::vector<uint64_t> get_node_colors(size_type v) const;  // 32 to 64
    bool is_colored(uint64_t);
    bool is_contained(dbg_boss::size_type v, uint64_t color) const;
    uint64_t get_color_num(size_type v) const;
    uint64_t get_sibling_color_num(size_type v) const;
    uint64_t get_value_vec_size();
    uint64_t get_used_value_vec_size();
    double get_color_vec_space();
    double get_color_matrix_space();
    double get_index_space();
    double get_ave_m_by_n();
    // void print_max_rss(std::string s="");
    void print_bloom_filter_info();

    size_type rev_comp(size_type v) const;
    range_t get_edges(size_type v) const;
    label_t node_label(size_type v) const;
    std::string node2string(size_type v, bool rev=true) const;
    size_type string2node(std::string string) const;
    size_type string2node(const char * string, size_t size) const;
    size_type tot_colored_nodes() const;
    uint64_t tot_colors() const;
    uint64_t tot_sibling_colors() const;
    double color_space_frac()const;

    std::pair<std::map<uint64_t, uint64_t>,
              std::map<uint64_t, uint64_t>> color_stats() const;


    std::vector<bs_t> backward_search(label_t &query, uint8_t mismatches,
                                      std::vector<edges_t::value_type> &syms=cs,
                                      std::vector<edges_t::size_type> &r_b=c_i,
                                      std::vector<edges_t::size_type> &r_a=c_j) const;
    std::vector<bs_t> backward_search(std::string query, uint8_t mismatches,
                                      std::vector<edges_t::value_type> &syms=cs,
                                      std::vector<edges_t::size_type> &r_b=c_i,
                                      std::vector<edges_t::size_type> &r_a=c_j) const;
    size_type tot_edges() const;
    size_type tot_nodes() const;
    void swap_color_matrix(color_matrix_t& other);

    //storage functions
    size_type serialize(std::ostream& out, structure_tree_node* v, std::string name) const;
    void load(std::istream&);

};

inline dbg_boss::range_t dbg_boss::get_edges(size_type v) const {
    size_t start, end;
    if(v==0){
        start=0;
    }else {
        start = m_node_marks_ss.select(v) + 1;
    }
    end = m_node_marks_ss.select(v+1);
    return std::make_pair(start, end);
}

inline dbg_boss::size_type dbg_boss::incomming(dbg_boss::size_type v, size_t value, bool is_symbol) const{
    assert(value>0);
    if(v==0) return 0;

    auto inv_res = m_f_array.inv_select(v);
    uint8_t symbol = inv_res.first;
    size_type sym_rank = inv_res.second;

    uint8_t bwt_sym = symbol<<1U;
    uint8_t m_bwt_sym = bwt_sym|1U;
    size_type bwt_pos = m_edge_bwt.select(sym_rank, bwt_sym);
    size_type next_bwt_pos;

    if(sym_rank==m_f_array.symbol_freq(symbol)){
        next_bwt_pos = m_edge_bwt.size();
    }else{
        next_bwt_pos = m_edge_bwt.select(sym_rank+1, bwt_sym);
    }

    size_type first_incoming = m_node_marks_rs.rank(bwt_pos);
    size_type marked_before = m_edge_bwt.rank(bwt_pos, m_bwt_sym);
    size_type n_marked = m_edge_bwt.rank(next_bwt_pos, m_bwt_sym)-
                         marked_before;
    if(is_symbol){
        label_t node_lab = node_label(first_incoming);
        if(node_lab.back()==dna_alphabet::char2comp[value]){
            return first_incoming;
        }

        size_type next_incomming;
        for(size_t i=1;i<=n_marked;i++){
            next_incomming =  m_node_marks_rs.rank(m_edge_bwt.select(marked_before+i,
                                                                     m_bwt_sym));
            node_lab = node_label(next_incomming);
            if(node_lab.back()==dna_alphabet::char2comp[value]){
                return next_incomming;
            }
        }
    }else{
        if(value==1){
            return first_incoming;
        }else{
            assert((value-1)<=n_marked);
            return m_node_marks_rs.rank(m_edge_bwt.select(marked_before+value-1, m_bwt_sym));
        }
    }
    return 0;
}

inline dbg_boss::size_type dbg_boss::outgoing(dbg_boss::size_type v, uint8_t value, bool is_symbol,
                                              std::vector<edges_t::value_type> &syms,
                                              std::vector<edges_t::size_type> &r_b,
                                              std::vector<edges_t::size_type> &r_a) const {

    range_t edges = get_edges(v);
    if(v<m_f_array.C[2]) return 0;

    if(!is_symbol){//value is treated as rank within the edges

        assert(value>0 && value<=(edges.second-edges.first+1));
        auto inv_select = m_edge_bwt.inverse_select(edges.first+value-1);
        uint8_t symbol = inv_select.second;
        size_type rank = inv_select.first;

        if(!(symbol & 1UL)){ //unmarked symbol
            return m_f_array.C[(symbol>>1U)] + rank;
        }else{
            return m_f_array.C[(symbol>>1U)] +
                   m_edge_bwt.rank(edges.first+value-1, (symbol & bit_clear)) - 1;
        }

    }else{//value is treated as a symbol
        value = dna_alphabet::char2comp[value];
        size_type n_symbols=0;
        m_edge_bwt.interval_symbols(edges.first, edges.second+1, n_symbols, syms, r_b, r_a);

        for(size_t i=0;i<n_symbols;i++){
            if((syms[i]>>1U) == value){
                if(!(syms[i] & 1U)){
                    return m_f_array.C[syms[i]>>1U] + r_b[i];
                }else{
                    return m_f_array.C[syms[i]>>1U] +
                           m_edge_bwt.rank(edges.first, (syms[i] & bit_clear)) - 1;
                }
            }
        }
        return 0;
    }
}//t is the rank inside the range [1..\sigma]

inline dbg_boss::size_type dbg_boss::tot_edges() const {
    return m_edge_bwt.size();
}

inline dbg_boss::size_type dbg_boss::tot_nodes() const {
    return m_f_array.C[6];
}

inline size_t dbg_boss::outdegree(dbg_boss::size_type v) const {
    range_t range = get_edges(v);
    return range.second-range.first+1;
}

inline size_t dbg_boss::indegree(dbg_boss::size_type v) const {

    uint8_t symbol = m_f_array[v];
    size_type symbol_freq = m_f_array.symbol_freq(symbol);
    size_type symbol_rank = v - m_f_array.C[symbol]+1;

    symbol = symbol<<1U;
    uint8_t symbol_marked = symbol | 1U;

    size_type edge_pos = m_edge_bwt.select(symbol_rank, symbol);
    size_type next_edge_pos;

    if(symbol_rank+1>symbol_freq){
        next_edge_pos = m_edge_bwt.size();
    }else{
        next_edge_pos = m_edge_bwt.select(symbol_rank+1, symbol);
    }

    size_type marked_bf = m_edge_bwt.rank(edge_pos, symbol_marked);
    size_type n_marked = m_edge_bwt.rank(next_edge_pos, symbol_marked)-marked_bf;

    return 1+n_marked;
}

inline std::vector<uint64_t> dbg_boss::get_node_colors(dbg_boss::size_type v) const {  // 32 to 64
    return m_color_matrix.get_colors(v);
}

inline bool dbg_boss::is_colored(uint64_t v) {
    return m_color_matrix.is_colored(v);
}

inline bool dbg_boss::is_contained(dbg_boss::size_type v, uint64_t color) const {
    return m_color_matrix.is_contained(v, color);
}

inline uint64_t dbg_boss::get_color_num(dbg_boss::size_type v) const {
    return m_color_matrix.get_color_num(v);
}

inline uint64_t dbg_boss::get_sibling_color_num(dbg_boss::size_type v) const {
    if (!m_color_matrix.is_colored(v)) return 0;

    size_type predecessor, sibling;
    size_t ind, outd;
    uint64_t color_num = 0;
    std::map<size_type, bool> siblings;

    ind = indegree(v);
    for (size_t i = 1; i <= ind; ++i) {
        predecessor = incomming(v, i, false);
        outd = outdegree(predecessor);
        for (size_t j = 1; j <= outd; ++j) {
            sibling = outgoing(predecessor, j);
            if (sibling == v) continue;
            siblings[sibling] = true;
        }
    }

    for (auto const &sibl : siblings) {
        color_num += get_color_num(sibl.first);
    }

    return color_num;
}

inline uint64_t dbg_boss::get_value_vec_size() {
    return m_color_matrix.get_value_vec_size();
}

inline uint64_t dbg_boss::get_used_value_vec_size() {
    return m_color_matrix.get_used_value_vec_size();
}

inline double dbg_boss::get_color_vec_space() {
    return m_color_matrix.get_color_vec_space();
}

inline double dbg_boss::get_color_matrix_space() {
    return size_in_mega_bytes(m_color_matrix);
}

inline double dbg_boss::get_index_space() {
    return size_in_mega_bytes(*this);
}

inline double dbg_boss::get_ave_m_by_n() {
    return m_color_matrix.get_ave_m_by_n();
}

// inline void dbg_boss::print_max_rss(std::string s) {
//     // getrusage関数による最大消費メモリ量のチェック
//     int chk;
//     struct rusage usage;
//     chk = getrusage(RUSAGE_SELF, &usage);
//     if(chk != 0){
//         printf("error\n");
//         exit(-1);
//     }
//     // このプロセスが消費した最大メモリ領域
//     printf("%sMax RSS = %lf MB\n", s.c_str(), usage.ru_maxrss / 1024.0);
// }

inline void dbg_boss::print_bloom_filter_info() {
    m_color_matrix.print_bloom_filter_info();
}

inline dbg_boss::size_type dbg_boss::tot_colored_nodes() const {
    return m_color_matrix.tot_colored_rows();
}

inline uint64_t dbg_boss::tot_colors() const {
    return m_color_matrix.tot_colors();
}

inline uint64_t dbg_boss::tot_sibling_colors() const {
    uint64_t count = 0, nei_col_count;

    for (size_type tmp_node = 0; tmp_node < tot_nodes(); ++tmp_node) {
        if (!m_color_matrix.is_colored(tmp_node)) continue;
        if (dna_alphabet::comp2char[node_label(tmp_node)[m_k - 2]] == '$') nei_col_count = 0;
        else nei_col_count = get_sibling_color_num(tmp_node);
        count += nei_col_count;
    }

    return count;
}

inline double dbg_boss::color_space_frac() const {
    return size_in_mega_bytes(m_color_matrix)/size_in_mega_bytes(*this);
}

inline std::pair<std::map<uint64_t, uint64_t>,
                 std::map<uint64_t, uint64_t>> dbg_boss::color_stats() const {
    return m_color_matrix.matrix_stats();
}

inline void dbg_boss::swap_color_matrix(color_matrix_t& other_cm) {
    m_color_matrix.swap(other_cm);
}

#endif //COL_BOSS_HPP
