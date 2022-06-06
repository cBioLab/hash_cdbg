#ifndef COL_BOSS_COLOR_MATRIX_HPP
#define COL_BOSS_COLOR_MATRIX_HPP


#include <sdsl/rrr_vector.hpp>
#include "bloom_filter.hpp"

template <class arr_t,
          class bv_t,
          class boss_t>
class color_matrix {

private:
    typedef sdsl::rrr_vector<63> comp_bv_t;
    typedef typename sdsl::rrr_vector<63>::rank_1_type comp_bv_rs_t;
    typedef typename sdsl::rrr_vector<63>::select_1_type comp_bv_ss_t;

    comp_bv_t m_row_marks;//rows are nodes
    comp_bv_rs_t m_row_marks_rs;
    comp_bv_ss_t m_row_marks_ss;
    comp_bv_t m_col_marks;//columns are the range of colors for a node
    comp_bv_ss_t m_col_marks_ss;
    bloom_filter<arr_t, boss_t> m_values;

public:

    typedef std::map<size_t, size_t> stat_map_t;
    typedef uint64_t size_type;

    struct matrix_skeleton{
        comp_bv_t colored_rows;
        comp_bv_t::rank_1_type colored_rows_rs;
        comp_bv_t colored_cells;
        comp_bv_t::select_1_type colored_cells_ss;
        sdsl::int_vector<64> cell_values;  // 16 to 64
        // std::map<uint64_t, std::map<uint64_t, bool> > neighbor_nodes;
    };

    bool is_colored(uint64_t v) const;
    std::vector<uint64_t> get_colors(uint64_t v) const;  // 32 to 64
    bool is_contained(uint64_t v, uint64_t color) const;
    uint64_t get_color_num(uint64_t v) const;
    uint64_t get_value_vec_size() const;
    uint64_t get_used_value_vec_size() const;
    double get_color_vec_space();
    uint64_t tot_colors() const;
    double get_ave_m_by_n() const;
    void print_bloom_filter_info() const;
    color_matrix()=default;
    color_matrix(arr_t& input_array,
                 bv_t& row_marks,
                 bv_t& col_marks);
    // explicit color_matrix(matrix_skeleton& m_esq, size_t n_threads);
    explicit color_matrix(matrix_skeleton& m_esq, boss_t &dbg_index, size_t n_threads);
    uint64_t tot_colored_rows() const;
    uint64_t tot_columns() const;  // 32 to 64
    std::pair<std::map<uint64_t, uint64_t>,
              std::map<uint64_t, uint64_t>> matrix_stats()const;

    uint64_t serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
    void load(std::istream&);
    void swap(color_matrix<arr_t, bv_t, boss_t>& other);
};

template<class arr_t, class bv_t, class boss_t>
inline bool color_matrix<arr_t, bv_t, boss_t>::is_colored(uint64_t v) const {
    return m_row_marks[v];
}

template<class arr_t, class bv_t, class boss_t>
color_matrix<arr_t, bv_t, boss_t>::color_matrix(arr_t&input_array,
                                        bv_t& row_marks,
                                        bv_t& col_marks) {

    comp_bv_t tmp_r_bv(row_marks);
    comp_bv_t tmp_c_bv(col_marks);

    m_row_marks.swap(tmp_r_bv);
    m_row_marks_rs.set_vector(&m_row_marks);
    m_col_marks.swap(tmp_c_bv);
    m_col_marks_ss.set_vector(&m_col_marks);

    bloom_filter<arr_t, boss_t> tmp_ef(input_array);
    m_values.swap(tmp_ef);
}

template<class arr_t, class bv_t, class boss_t>
inline std::vector<uint64_t> color_matrix<arr_t, bv_t, boss_t>::get_colors(uint64_t v) const{  // 32 to 64

    std::vector<uint64_t> colors;  // 32 to 64

    if(!m_row_marks[v]) return colors;

    uint64_t row = m_row_marks_rs(v);
    uint64_t col_start = m_col_marks_ss(row+1);
    uint64_t col_end = m_col_marks_ss(row+2)-1;
    uint64_t acc = 0;  // size_t to uint64_t

    for(uint64_t i=col_start; i<=col_end;i++){
        acc+=m_values.read_elm(i);
        colors.push_back(acc);
    }

    return colors;
}

template<class arr_t, class bv_t, class boss_t>
inline bool color_matrix<arr_t, bv_t, boss_t>::is_contained(uint64_t v, uint64_t color) const {
    if (!m_row_marks[v]) return false;  // 色を付けるノードじゃなかった場合

    uint64_t row = m_row_marks_rs(v);  // ノードが何番目の色付きノードかを取得 (0~)
    uint64_t col_num = m_col_marks_ss(row + 2) - m_col_marks_ss(row + 1);

    return m_values.is_contained(row, col_num, color);
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::get_color_num(uint64_t v) const {
    if (!m_row_marks[v]) return 0;  // 色を付けるノードじゃなかった場合

    uint64_t row = m_row_marks_rs(v);
    uint64_t col_start = m_col_marks_ss(row + 1);
    uint64_t col_end = m_col_marks_ss(row + 2) - 1;

    return col_end - col_start + 1;
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::get_value_vec_size() const {
    return m_values.get_value_vec_size();
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::get_used_value_vec_size() const {
    return m_values.get_used_value_vec_size();
}

template<class arr_t, class bv_t, class boss_t>
inline double color_matrix<arr_t, bv_t, boss_t>::get_color_vec_space() {
    return size_in_mega_bytes(m_values);
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::tot_colors() const {
    return m_col_marks_ss(tot_colored_rows() + 1);
}

template<class arr_t, class bv_t, class boss_t>
inline double color_matrix<arr_t, bv_t, boss_t>::get_ave_m_by_n() const {
    uint64_t len = get_used_value_vec_size();
    uint64_t count = tot_colors();
    return (double)len / (double)count;
}

template<class arr_t, class bv_t, class boss_t>
inline void color_matrix<arr_t, bv_t, boss_t>::print_bloom_filter_info() const {
    m_values.print_bloom_filter_info();
}

template<class arr_t, class bv_t, class boss_t>
void color_matrix<arr_t, bv_t, boss_t>::swap(color_matrix<arr_t, bv_t, boss_t>& other) {

    m_row_marks.swap(other.m_row_marks);
    m_row_marks_rs.swap(other.m_row_marks_rs);
    m_row_marks_rs.set_vector(&m_row_marks);

    m_col_marks.swap(other.m_col_marks);
    m_col_marks_ss.swap(other.m_col_marks_ss);
    m_col_marks_ss.set_vector(&m_col_marks);

    m_values.swap(other.m_values);
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::tot_colored_rows() const {
    return m_row_marks_rs(m_row_marks.size());
}

template<class arr_t, class bv_t, class boss_t>
inline uint64_t color_matrix<arr_t, bv_t, boss_t>::tot_columns() const {  // 32 to 64
    return m_values.max_val();
}

template<class arr_t, class bv_t, class boss_t>
uint64_t color_matrix<arr_t, bv_t, boss_t>::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;

    written_bytes += m_row_marks.serialize(out,child, "m_row_marks");
    written_bytes += m_row_marks_rs.serialize(out, child, "m_row_marks_rs");
    written_bytes += m_col_marks.serialize(out, child, "m_col_marks");
    written_bytes += m_col_marks_ss.serialize(out, child, "m_col_marks_ss");
    written_bytes += m_values.serialize(out, child, "m_values");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class arr_t, class bv_t, class boss_t>
void color_matrix<arr_t, bv_t, boss_t>::load(std::istream &in) {
    m_row_marks.load(in);
    m_row_marks_rs.load(in, &m_row_marks);
    m_col_marks.load(in);
    m_col_marks_ss.load(in, &m_col_marks);
    m_values.load(in);
}

template<class arr_t, class bv_t, class boss_t>
std::pair<std::map<uint64_t, uint64_t>,
          std::map<uint64_t, uint64_t>> color_matrix<arr_t, bv_t, boss_t>::matrix_stats() const {

    std::map<uint64_t, uint64_t> colored_node_freqs;
    std::map<uint64_t, uint64_t> color_freqs;
    size_type start, end;
    size_type tot_col_rows = tot_colored_rows();

    for(size_type i=1;i<=tot_col_rows;i++){

        start = m_col_marks_ss(i);
        end = m_col_marks_ss(i+1)-1;
        colored_node_freqs[end-start+1]++;

        for(size_type j=start;j<=end;j++){
            color_freqs[m_values.read_elm(j)]++;
        }
    }
    return {color_freqs, colored_node_freqs};
}

template<class arr_t, class bv_t, class boss_t>
// color_matrix<arr_t, bv_t, boss_t>::color_matrix(color_matrix::matrix_skeleton &m_esq, size_t n_threads) {
color_matrix<arr_t, bv_t, boss_t>::color_matrix(color_matrix::matrix_skeleton &m_esq, boss_t &dbg_index, size_t n_threads) {
    m_row_marks.swap(m_esq.colored_rows);
    m_row_marks_rs.set_vector(&m_row_marks);
    m_row_marks_ss.set_vector(&m_row_marks);

    m_col_marks.swap(m_esq.colored_cells);
    m_col_marks_ss.set_vector(&m_col_marks);

    // bloom_filter<arr_t, boss_t> tmp_ef(m_esq.cell_values);
    // bloom_filter<arr_t, boss_t> tmp_ef(m_esq.cell_values, m_esq.neighbor_nodes, m_col_marks_ss, tot_colored_rows(), n_threads);
    bloom_filter<arr_t, boss_t> tmp_ef(m_esq.cell_values, std::ref(dbg_index), m_row_marks_rs, m_row_marks_ss, m_col_marks_ss, tot_colored_rows(), n_threads);

    m_values.swap(tmp_ef);
}

#endif //COL_BOSS_COLOR_MATRIX_HPP
