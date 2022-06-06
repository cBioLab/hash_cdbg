#ifndef COL_BOSS_BUILD_INDEX_HPP
#define COL_BOSS_BUILD_INDEX_HPP

#include <iostream>
#include <thread>
#include <mutex>
#include <stdexcept>
#include "f_array_t.hpp"
#include "fastx_parser.hpp"
#include "dna_alphabet.hpp"
#include "boss.hpp"
#include "color_matrix.hpp"


class build_index {
private:
    typedef int_vector_buffer<8>              text_t;
    typedef dbg_boss::size_type               size_type;
    typedef dbg_boss::color_matrix_t          color_matrix_t;
    typedef color_matrix_t::matrix_skeleton   matrix_skeleton;

public:
    static void build_boss(std::string input_file, sdsl::cache_config &config, size_t K);
    static color_matrix_t color_dbg(dbg_boss &dbg_index, sdsl::cache_config &config, size_t n_theads);

private:

    static void color_dbg_internal(dbg_boss &dbg_index,
                                   long long int start,
                                   long long int end,
                                   matrix_skeleton& m_skeleton,
                                   sdsl::cache_config &config);

    static void build_matrix_skeleton(dbg_boss &dbg_index,
                                      matrix_skeleton &skeleton,
                                      size_t n_threads,
                                      sdsl::cache_config &config);

    static void estimate_color_rows(dbg_boss &dbg_index,
                                    size_type start,
                                    size_type end,
                                    sdsl::bit_vector &color_marks);

    static void estimate_color_cells(dbg_boss &dbg_index,
                                     long long int start,
                                     long long int end,
                                     sdsl::cache_config &config,
                                     sdsl::rrr_vector<63> &color_marks,
                                     sdsl::rrr_vector<63>::rank_1_type &color_marks_rs,
                                    //  std::map<build_index::size_type, std::map<build_index::size_type, bool> > &n_nodes,
                                     std::map<size_type, uint64_t> &c_nodes);  // 32 to 64

    static void estimate_amb_seqs(dbg_boss &dbg_index,
                                  long long int start,
                                  long long int end,
                                  sdsl::cache_config &config,
                                  sdsl::rrr_vector<63> &color_marks,
                                  sdsl::rrr_vector<63>::rank_1_type &color_marks_rs,
                                  size_t& tot_amb);


    static void build_SA_BWT_LCP(sdsl::cache_config &config);
    static void build_edge_BWT(sdsl::cache_config &config, size_t K);
    static void build_KLCP(sdsl::cache_config &config, size_t K);

};


#endif //COL_BOSS_BUILD_INDEX_HPP
