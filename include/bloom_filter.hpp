#ifndef COL_BOSS_BLOOM_FILTER_HPP
#define COL_BOSS_BLOOM_FILTER_HPP

#include <sdsl/rrr_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <algorithm>
#include <bitset>
#include <thread>
#include <unordered_set>

template<class arr_t, class boss_t>
class bloom_filter {

public:
    typedef uint64_t size_type;

private:
    typedef sdsl::rrr_vector<63> comp_bv_t;
    typedef typename sdsl::rrr_vector<63>::rank_1_type comp_bv_rs_t;
    typedef typename sdsl::rrr_vector<63>::select_1_type comp_bv_ss_t;

    comp_bv_t m_values_marks;
    comp_bv_ss_t m_values_marks_ss;
    sdsl::int_vector<64> m_values;
    size_t m_w=64;
    uint64_t m_size=0;
    uint64_t m_max_val=0;  // 32 to 64
    constexpr static double fp_rate = 0.00001;

public:
    explicit bloom_filter(arr_t& uncomp_values);
    // explicit bloom_filter(arr_t &uncomp_values, std::map<uint64_t, std::map<uint64_t, bool> > &neighbor_nodes, comp_bv_ss_t &col_marks_ss, uint64_t tot_colored_rows, size_t n_threads);
    explicit bloom_filter(arr_t &uncomp_values, boss_t &dbg_index, comp_bv_rs_t &row_marks_rs, comp_bv_ss_t &row_marks_ss, comp_bv_ss_t &col_marks_ss, uint64_t tot_colored_rows, size_t n_threads);
    bloom_filter()=default;
    size_type size() const;
    size_type max_val() const;
    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const;
    void load(std::istream&);
    size_type read_elm(size_type pos) const;
    void swap(bloom_filter<arr_t, boss_t>& other);
    static uint64_t calc_max_prime(uint64_t n);
    bool is_contained(uint64_t row, uint64_t col_num, uint64_t value) const;
    uint64_t get_value_vec_size() const;
    uint64_t get_used_value_vec_size() const;
    void print_bloom_filter_info() const;

private:
    void bits_write(uint64_t i, uint64_t j, uint64_t value, sdsl::int_vector_buffer<64>& buffer);
    size_type bits_read(size_type i, size_type j) const;
    static void build_filter(uint64_t hash_num, uint64_t filter_size, uint64_t value, sdsl::int_vector<1> &filter);
    // static void build_filters(arr_t &uncomp_values, std::map<uint64_t, std::map<uint64_t, bool> > &neighbor_nodes, comp_bv_ss_t &col_marks_ss, size_t start, size_t end, std::map<uint64_t, uint64_t> &m_by_n_map, std::map<uint64_t, sdsl::int_vector<1> > &filter_map);
    // static void build_filters(arr_t &uncomp_values, boss_t &dbg_index, comp_bv_rs_t &row_marks_rs, comp_bv_ss_t &row_marks_ss, comp_bv_ss_t &col_marks_ss, size_t start, size_t end, std::map<uint64_t, sdsl::int_vector<1> > &filter_map, std::map<uint64_t, uint64_t> &sibl_col_rate_map);
    static void build_filters(arr_t &uncomp_values, boss_t &dbg_index, comp_bv_rs_t &row_marks_rs, comp_bv_ss_t &row_marks_ss, comp_bv_ss_t &col_marks_ss, size_t start, size_t end, std::map<uint64_t, uint64_t> &col_num_map, std::map<uint64_t, sdsl::int_vector<1> > &filter_map, std::map<uint64_t, uint64_t> &sibl_col_rate_map);
    static uint64_t hash1(uint64_t value, uint64_t filter_size);
    static uint64_t hash2(uint64_t value, uint64_t filter_size);
    static uint64_t hash_bloom(uint64_t h1, uint64_t h2, uint64_t filter_size, uint64_t i);
    static uint64_t calc_hash_num(uint64_t filter_size, uint64_t col_num);
    void set_filter(uint64_t start, uint64_t filter_size, sdsl::int_vector<1> filter, sdsl::int_vector_buffer<64>& values_buffer);
    static bool check_contained(uint64_t hash_num, uint64_t filter_size, uint64_t value, sdsl::int_vector<1> &filter);
};

template<class arr_t, class boss_t>
bloom_filter<arr_t, boss_t>::bloom_filter(arr_t &uncomp_values) {

    size_t width;
    uint64_t curr_pos=0, value;
    m_size = uncomp_values.size();

    sdsl::int_vector_buffer<64> values_buffer("ef_buffer.sdsl",
                                              std::ios::out,
                                              1024*1024,
                                              false);

    sdsl::int_vector_buffer<1> bv_buffer("ef_buffer_bv.sdsl",
                                         std::ios::out,
                                         1024*1024,
                                         false);

    for(size_type i=0;i<uncomp_values.size();i++){

        //2 is the smallest number allowed
        value = uncomp_values[i]+2;

        if((value-2)>m_max_val) m_max_val = value-2;

        width =  (m_w-1)-__builtin_clzll(value);
        bits_write(curr_pos,
                   curr_pos+width-1,
                   value & ~(1ULL << width),
                   values_buffer);
        bv_buffer[curr_pos] = true;
        curr_pos+=width;
    }
    bv_buffer[curr_pos] = true;

    /*std::cout<<"values"<<std::endl;
    for(size_type i=0;i<values_buffer.size();i++){
        std::cout<<values_buffer[i]<<" ";
    }
    std::cout<<" "<<std::endl;*/

    values_buffer.close();
    bv_buffer.close();

    sdsl::int_vector<64> values;
    sdsl::load_from_file(m_values, "ef_buffer.sdsl");

    sdsl::bit_vector tmp_bv;
    sdsl::load_from_file(tmp_bv, "ef_buffer_bv.sdsl");

    sdsl::rrr_vector<63> tmp_bv_comp(tmp_bv);
    m_values_marks.swap(tmp_bv_comp);

    /*for(size_type i=0;i<m_values.size();i++){
        std::cout<<m_values[i];
    }
    std::cout<<""<<std::endl;*/

    m_values_marks_ss.set_vector(&m_values_marks);

    if(remove("ef_buffer.sdsl") != 0) {
        perror("Error deleting file");
    }

    if(remove("ef_buffer_bv.sdsl") != 0) {
        perror("Error deleting file");
    }
}

template<class arr_t, class boss_t>
inline bool bloom_filter<arr_t, boss_t>::check_contained(uint64_t hash_num,
                                             uint64_t filter_size,
                                             uint64_t value,
                                             sdsl::int_vector<1> &filter) {

    uint64_t pos;
    uint64_t h1 = hash1(value, filter_size);
    uint64_t h2 = hash2(value, filter_size);
    for (size_type i = 0; i < hash_num; ++i) {
        pos = hash_bloom(h1, h2, filter_size, i);
        if (!filter[pos]) return false;
    }

    return true;
}

template<class arr_t, class boss_t>
inline void bloom_filter<arr_t, boss_t>::build_filter(uint64_t hash_num,
                                        uint64_t filter_size,
                                        uint64_t value,
                                        sdsl::int_vector<1> &filter) {

    uint64_t pos;
    uint64_t h1 = hash1(value, filter_size);
    uint64_t h2 = hash2(value, filter_size);
    for (size_type i = 0; i < hash_num; ++i) {
        pos = hash_bloom(h1, h2, filter_size, i);
        filter[pos] = 1;
    }
}

template<class arr_t, class boss_t>
void bloom_filter<arr_t, boss_t>::build_filters(arr_t &uncomp_values,
                                    // std::map<uint64_t, std::map<uint64_t, bool> > &neighbor_nodes,
                                    boss_t &dbg_index,
                                    comp_bv_rs_t &row_marks_rs,
                                    comp_bv_ss_t &row_marks_ss,
                                    comp_bv_ss_t &col_marks_ss,
                                    size_t start,
                                    size_t end,
                                    std::map<uint64_t, uint64_t> &col_num_map,
                                    std::map<uint64_t, sdsl::int_vector<1> > &filter_map,
                                    std::map<uint64_t, uint64_t> &sibl_col_rate_map) {

    uint64_t col_start, col_end, col_num;
    uint64_t nei_start, nei_end;
    uint64_t filter_size, hash_num, all_count, fp_count;
    uint64_t filter_size_start, filter_size_add;
    uint64_t limit;
    size_t ind, outd;
    size_type prev_node, nei_node;

    for (size_t i = start; i < end; ++i) {
        col_start = col_marks_ss(i + 1);
        col_end = col_marks_ss(i + 2) - 1;
        col_num = col_end - col_start + 1;
        filter_size_start = 2;
        filter_size_add = 1;
        limit = col_num * 300;

        std::vector<std::pair<uint64_t, uint64_t> > sibling_color_idxs;
        uint64_t sibling_color_num = 0;
        if (dbg_index.node2string(row_marks_ss(i + 1))[0] != '$') {
            ind = dbg_index.indegree(row_marks_ss(i + 1));
            for (size_t j = 1; j <= ind; ++j) {
                prev_node = dbg_index.incomming(row_marks_ss(i + 1), j, false);
                outd = dbg_index.outdegree(prev_node);
                if(outd > 1) {
                    for (size_t k = 1; k <= outd; ++k) {
                        nei_node = dbg_index.outgoing(prev_node, k, false);
                        if (nei_node == row_marks_ss(i + 1)) continue;
                        nei_start = col_marks_ss(row_marks_rs(nei_node) + 1);
                        nei_end = col_marks_ss(row_marks_rs(nei_node) + 2) - 1;
                        sibling_color_idxs.emplace_back(std::make_pair(nei_start, nei_end));
                        sibling_color_num += nei_end - nei_start + 1;
                    }
                }
            }
        }
        sibl_col_rate_map[i] = sibling_color_num / col_num;
        if (sibling_color_num == 0) {
            sdsl::bit_vector tmp_filter;
            sdsl::util::assign(tmp_filter, sdsl::int_vector<1>(1, 1));
            filter_map[i] = tmp_filter;
            col_num_map[i] = col_num;
            continue;
        }
        // if (sibling_color_num <= col_num * 2) {
        //     filter_size_start = 2;
        //     filter_size_add = 1;
        // } else {
        //     filter_size_start = col_num;
        //     filter_size_add = col_num;
        // }

        std::unordered_set<uint64_t> own_colors;
        for (size_type j = col_start; j <= col_end; ++j) {
            own_colors.insert(uncomp_values[j]);
        }

        for (filter_size = filter_size_start; filter_size < limit; filter_size += filter_size_add) {
            hash_num = calc_hash_num(filter_size, col_num);
            sdsl::bit_vector tmp_filter;
            sdsl::util::assign(tmp_filter, sdsl::int_vector<1>(filter_size, 0));
            // for (size_type j = col_start; j <= col_end; ++j) {
            //     build_filter(hash_num, filter_size, uncomp_values[j], std::ref(tmp_filter));
            // }
            for (auto const &color : own_colors) {
                build_filter(hash_num, filter_size, color, std::ref(tmp_filter));
            }

            all_count = 0;
            fp_count = 0;
            // for (auto const &nei : neighbor_nodes[i]) {
            //     nei_start = col_marks_ss(nei.first + 1);
            //     nei_end = col_marks_ss(nei.first + 2) - 1;
            //     for (size_type k = nei_start; k <= nei_end; ++k) {
            //         all_count++;
            //         if (check_contained(hash_num, filter_size, uncomp_values[k], std::ref(tmp_filter))) {
            //             fp_count++;
            //         }
            //     }
            // }

            for (auto const &[sibl_col_start, sibl_col_end] : sibling_color_idxs) {
                for (size_type j = sibl_col_start; j <= sibl_col_end; ++j) {
                    // bool corrision_flag = false;
                    // for (size_type k = col_start; k <= col_end; ++k) {
                    //     if (uncomp_values[j] == uncomp_values[k]) {
                    //         corrision_flag = true;
                    //         break;
                    //     }
                    // }
                    // if (corrision_flag) continue;
                    if (own_colors.find(uncomp_values[j]) != own_colors.end()) continue;
                    all_count++;
                    if (check_contained(hash_num, filter_size, uncomp_values[j], std::ref(tmp_filter))) {
                        fp_count++;
                    }
                }
            }

            if (all_count == 0 || (double)fp_count / (double)all_count <= fp_rate || filter_size + filter_size_add >= limit) {
                filter_map[i] = tmp_filter;
                col_num_map[i] = col_num;
                own_colors.clear();
                // sibling_colors.clear();
                break;
            }
            // if (filter_size / col_num >= 2) filter_size_add = col_num;
        }
    }
}

template<class arr_t, class boss_t>
bloom_filter<arr_t, boss_t>::bloom_filter(arr_t &uncomp_values,
                        //   std::map<uint64_t, std::map<uint64_t, bool> > &neighbor_nodes,
                          boss_t &dbg_index,
                          comp_bv_rs_t &row_marks_rs,
                          comp_bv_ss_t &row_marks_ss,
                          comp_bv_ss_t &col_marks_ss,
                          uint64_t tot_colored_rows,
                          size_t n_threads) {

    uint64_t curr_pos = 0, filter_size;
    m_size = uncomp_values.size();

    sdsl::int_vector_buffer<64> values_buffer("ef_buffer.sdsl",
                                              std::ios::out,
                                              1024*1024,
                                              false);

    sdsl::int_vector_buffer<1> bv_buffer("ef_buffer_bv.sdsl",
                                         std::ios::out,
                                         1024*1024,
                                         false);

    size_t start, end;
    std::vector<std::map<uint64_t, uint64_t> > col_num_maps(n_threads);
    std::vector<std::map<uint64_t, sdsl::int_vector<1> > > filter_maps(n_threads);
    std::vector<std::map<uint64_t, uint64_t> > sibl_col_rate_maps(n_threads);
    size_t elm_per_thread = tot_colored_rows / n_threads;
    size_t rem = tot_colored_rows % n_threads;
    std::vector<std::thread> filter_threads;
    std::cout << "fp rate: " << fp_rate << std::endl;  // for debug
    for (size_t i = 0; i < n_threads; ++i) {
        start = i * elm_per_thread;
        end = (i + 1) * elm_per_thread;
        if (i == n_threads - 1) end += rem;
        filter_threads.emplace_back(std::thread(build_filters,
                                                std::ref(uncomp_values),
                                                // std::ref(neighbor_nodes),
                                                std::ref(dbg_index),
                                                std::ref(row_marks_rs),
                                                std::ref(row_marks_ss),
                                                std::ref(col_marks_ss),
                                                start,
                                                end,
                                                std::ref(col_num_maps[i]),
                                                std::ref(filter_maps[i]),
                                                std::ref(sibl_col_rate_maps[i])));
    }
    for (auto &it : filter_threads) it.join();

    std::map<uint64_t, uint64_t> mbyn_statistic;
    // std::map<uint64_t, std::map<uint64_t, uint64_t> > mbyn_statistic_detail;
    std::map<uint64_t, std::map<uint64_t, uint64_t> > mbyn_siblcolrate_statistic;
    for (size_t i = 0; i < n_threads; ++i) {
        for (auto const &node : filter_maps[i]) {
            filter_size = node.second.size();
            set_filter(curr_pos, filter_size, node.second, std::ref(values_buffer));

            uint64_t col_num = col_num_maps[i][node.first];
            // hash_num = calc_hash_num(filter_size, col_num);
            bv_buffer[curr_pos] = true;
            // bv_buffer[curr_pos + hash_num] = true;
            curr_pos += filter_size;

            mbyn_statistic[(uint64_t)(filter_size / col_num)]++;
            // mbyn_statistic_detail[(uint64_t)(filter_size / col_num)][(uint64_t)(filter_size * 10 / col_num) % 10]++;
            mbyn_siblcolrate_statistic[(uint64_t)(filter_size / col_num)][sibl_col_rate_maps[i][node.first]]++;
        }
    }
    bv_buffer[curr_pos] = true;

    for (auto const &mbn : mbyn_statistic) {
        std::cout << mbn.first << ": " << mbn.second << std::endl;
        // for (auto const &mbn_detail : mbyn_statistic_detail[mbn.first]) {
        //     std::cout << "  0." << mbn_detail.first << ": " << mbn_detail.second << std::endl;
        // }
        for (auto const &tmp_col_rate : mbyn_siblcolrate_statistic[mbn.first]) {
            std::cout << "  " << tmp_col_rate.first << ": " << tmp_col_rate.second << std::endl;
        }
    }
    mbyn_statistic.clear();
    // mbyn_statistic_detail.clear();
    mbyn_siblcolrate_statistic.clear();

    // std::cout << "value_buffer" << std::endl;
    // for (size_type i = 0; i < values_buffer.size(); ++i) {
    //     std::bitset<64> bit(values_buffer[i]);
    //     for (size_type j = 0; j < 64; ++j) {
    //         std::cout << bit[j];
    //     }
    //     std::cout << " ";
    // }
    // std::cout << std::endl;

    // for (size_type i = 0; i < bv_buffer.size(); ++i) {
    //     std::cout << bv_buffer[i];
    //     if (i % 64 == 63) std::cout << " ";
    // }
    // std::cout << std::endl;

    values_buffer.close();
    bv_buffer.close();

    sdsl::int_vector<64> values;
    sdsl::load_from_file(m_values, "ef_buffer.sdsl");

    sdsl::bit_vector tmp_bv;
    sdsl::load_from_file(tmp_bv, "ef_buffer_bv.sdsl");

    sdsl::rrr_vector<63> tmp_bv_comp(tmp_bv);
    m_values_marks.swap(tmp_bv_comp);

    m_values_marks_ss.set_vector(&m_values_marks);

    if(remove("ef_buffer.sdsl") != 0) {
        perror("Error deleting file");
    }

    if(remove("ef_buffer_bv.sdsl") != 0) {
        perror("Error deleting file");
    }
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::calc_max_prime(uint64_t n) {
    std::vector<int> prime = { 1,  1,  1,  2,  3,  3,  5,  5,  7,  7,  7,  7, 11, 11, 13, 13, 13, 13, 17, 17,
                                  19, 19, 19, 19, 23, 23, 23, 23, 23, 23, 29, 29, 31, 31, 31, 31, 31, 31, 37, 37,
                                  37, 37, 41, 41, 43, 43, 43, 43, 47, 47, 47, 47, 47, 47, 53, 53, 53, 53, 53, 53,
                                  59, 59, 61, 61, 61, 61, 61, 61, 67, 67, 67, 67, 71, 71, 73, 73, 73, 73, 73, 73,
                                  79, 79, 79, 79, 83, 83, 83, 83, 83, 83, 89, 89, 89, 89, 89, 89, 89, 89, 97, 97};

    if (n < 100) return prime[n];
    if (n < 300) return 97;
    if (n < 600) return 293;
    if (n < 1000) return 599;
    if (n < 2000) return 997;
    if (n < 4000) return 1999;
    if (n < 6000) return 3989;
    if (n < 10000) return 5987;
    return 9973;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::hash1(uint64_t value, uint64_t filter_size) {
    return value % filter_size;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::hash2(uint64_t value, uint64_t filter_size) {
    uint64_t prime = calc_max_prime(filter_size);
    return prime - (value % prime);
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::hash_bloom(uint64_t h1, uint64_t h2, uint64_t filter_size, uint64_t i) {
    return (h1 + i * h2) % filter_size;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::calc_hash_num(uint64_t filter_size, uint64_t col_num) {
    return (uint64_t)std::max((double)(filter_size / col_num) * 0.7, 1.0);
}

template<class arr_t, class boss_t>
inline void bloom_filter<arr_t, boss_t>::set_filter(uint64_t start,
                                        uint64_t filter_size,
                                        sdsl::int_vector<1> filter,
                                        sdsl::int_vector_buffer<64>& values_buffer) {
    size_t cell_start, cell_end;
    cell_start = start / m_w;
    cell_end = (start + filter_size - 1) / m_w;
    if (cell_end >= values_buffer.size()) {
        for (size_type i = values_buffer.size(); i <= cell_end; ++i) {
            values_buffer[i] = 0;
        }
    }

    size_t curr_cell = cell_start;
    size_t curr_pos = start % m_w;
    for (size_t i = 0; i < filter_size; ++i) {
        if (filter[i]) {
            values_buffer[curr_cell] = values_buffer[curr_cell] | (1ULL << curr_pos);
        }
        curr_pos++;
        if (curr_pos == m_w) {
            curr_pos = 0;
            curr_cell++;
        }
    }
}

template<class arr_t, class boss_t>
inline bool bloom_filter<arr_t, boss_t>::is_contained(uint64_t row, uint64_t col_num, uint64_t value) const {
    uint64_t start = m_values_marks_ss(row + 1);
    // uint64_t hash_num = m_values_marks_ss(row * 2 + 2) - start;
    uint64_t end = m_values_marks_ss(row + 2) - 1;
    uint64_t filter_size = end - start + 1;
    assert(filter_size > 0);
    uint64_t hash_num = calc_hash_num(filter_size, col_num);

    uint64_t pos;
    uint64_t h1 = hash1(value, filter_size);
    uint64_t h2 = hash2(value, filter_size);

    // std::cout << "  h1: " << h1 << ", h2: " << h2 << ", ";
    for (size_type i = 0; i < hash_num; ++i) {
        pos = hash_bloom(h1, h2, filter_size, i) + start;
        // std::cout << pos << "(" << pos - start << ") ";
        if ((m_values[pos / m_w] & (1ULL << (pos & (m_w - 1)))) == 0){
            // std::cout << std::endl;
            return false;
        }
    }
    // std::cout << std::endl;

    return true;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::get_value_vec_size() const {
    return (uint64_t)m_values.size() * m_w;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::get_used_value_vec_size() const {
    return (uint64_t)m_values_marks.size() - 1;
}

template<class arr_t, class boss_t>
inline void bloom_filter<arr_t, boss_t>::print_bloom_filter_info() const {
    std::cout << "m_size           : " << m_size << std::endl;
    std::cout << "m_max_val        : " << m_max_val << std::endl;
    std::cout << "m_values         : " << get_value_vec_size() << " bits" << std::endl;
    std::cout << "m_values         : " << size_in_mega_bytes(m_values) << " MB" << std::endl;
    std::cout << "m_values_marks   : " << m_values_marks.size() << " bits" << std::endl;
    std::cout << "m_values_marks   : " << size_in_mega_bytes(m_values_marks) << " MB" << std::endl;

    sdsl::rrr_vector<63>::rank_1_type m_values_marks_rs;
    m_values_marks_rs.set_vector(&m_values_marks);
    std::cout << "m_values_marks 1 : " << m_values_marks_rs(m_values_marks.size()) << " bits" << std::endl;
}

template<class arr_t, class boss_t>
inline void bloom_filter<arr_t, boss_t>::bits_write(uint64_t i,
                                        uint64_t j,
                                        uint64_t value,
                                        sdsl::int_vector_buffer<64>& buffer) {

    size_t cell_i, cell_j;
    cell_i = i/m_w;
    cell_j = j/m_w;

    if(cell_i>=buffer.size()) buffer[cell_i] = 0;
    if(cell_i==cell_j){
        buffer[cell_j] = buffer[cell_j] & ~(((1ULL<<(j-i+1)) -1) << (i & (m_w-1)));  // valueの値を格納するとこを0にセット
        buffer[cell_j] = buffer[cell_j] | value << (i & (m_w-1));  // valueの値を格納
        //std::cout<<buffer[cell_j]<<std::endl;
    }else{
        if(cell_j>=buffer.size()) buffer[cell_j] = 0;
        buffer[cell_i] = (buffer[cell_i] & ((1ULL<<(i & (m_w-1)))-1)) | (value << (i & (m_w-1)));  // valueの下m_w-(i%m_w)桁を格納
        buffer[cell_j] = (buffer[cell_j] & ~((1ULL<<((j+1) & (m_w-1)))-1)) |
                           (value >> (m_w-(i & (m_w-1))));  // valueの残りを格納
    }
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::bits_read(bloom_filter::size_type i, bloom_filter::size_type j) const {
    size_t cell_i, cell_j;
    cell_i = i/m_w;
    cell_j = j/m_w;

    if(cell_i == cell_j){
        return (m_values[cell_j] >> (i & (m_w-1))) & ((1ULL<<(j-i+1))-1);
    }else{
        return (m_values[cell_i] >> (i & (m_w-1) )  |
                (m_values[cell_j] & ((1ULL <<((j+1) & (m_w-1)))-1)) << (m_w - (i & (m_w-1))));
    }
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::read_elm(bloom_filter::size_type pos) const {
    assert(pos<m_size);
    size_type i= m_values_marks_ss(pos+1);
    size_type j= m_values_marks_ss(pos+2)-1;
    size_type elm = bits_read(i,j);

    return (1ULL<<(j-i+1ULL) | elm)-2;
}

template<class arr_t, class boss_t>
uint64_t bloom_filter<arr_t, boss_t>::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += m_values.serialize(out,child, "m_values");
    written_bytes += m_values_marks.serialize(out, child, "m_values_marks");
    written_bytes += m_values_marks_ss.serialize(out, child, "m_values_marks_ss");
    written_bytes += write_member(m_size, out, child, "m_size");
    written_bytes += write_member(m_max_val, out, child, "m_max_val");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template<class arr_t, class boss_t>
void bloom_filter<arr_t, boss_t>::load(std::istream & in) {
    m_values.load(in);
    m_values_marks.load(in);
    m_values_marks_ss.load(in, &m_values_marks);
    sdsl::read_member(m_size, in);
    sdsl::read_member(m_max_val, in);
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::size() const {
    return m_size;
}

template<class arr_t, class boss_t>
inline uint64_t bloom_filter<arr_t, boss_t>::max_val() const {
    return m_max_val;
}

template<class arr_t, class boss_t>
void bloom_filter<arr_t, boss_t>::swap(bloom_filter<arr_t, boss_t> &other) {
    m_values.swap(other.m_values);
    m_values_marks.swap(other.m_values_marks);
    m_values_marks_ss.swap(other.m_values_marks_ss);
    m_values_marks_ss.set_vector(&m_values_marks);
    std::swap(m_size,other.m_size);
    std::swap(m_max_val,other.m_max_val);
}

#endif //COL_BOSS_BLOOM_FILTER_HPP
