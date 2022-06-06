#include <boss.hpp>
#include <build_index.hpp>

#include <sys/stat.h>
#include <bits/stdc++.h>
#include <sys/types.h>

std::vector<dbg_boss::edges_t::value_type> dbg_boss::cs(2*dna_alphabet::sigma);
std::vector<dbg_boss::edges_t::size_type> dbg_boss::c_i(2*dna_alphabet::sigma);
std::vector<dbg_boss::edges_t::size_type> dbg_boss::c_j(2*dna_alphabet::sigma);

dbg_boss::dbg_boss(std::string &input_file, size_t K, size_t n_threads) : m_n_solid_nodes(0),
                                                                          m_k(K){

    cache_config config;

    build_index::build_boss(input_file, config, K);

    load_from_file(m_edge_bwt, cache_file_name("ebwt_wt_huff",config));
    load_from_file(m_f_array, cache_file_name("f_array",config));

    load_from_file(m_node_marks, cache_file_name("node_marks", config));
    m_node_marks_ss.set_vector(&m_node_marks);
    m_node_marks_rs.set_vector(&m_node_marks);

    load_from_file(m_solid_nodes, cache_file_name("solid_nodes", config));
    m_solid_nodes_rs.set_vector(&m_solid_nodes);
    m_solid_nodes_ss.set_vector(&m_solid_nodes);
    m_n_solid_nodes = m_solid_nodes_rs(m_solid_nodes.size());

    auto tmp_color_matrix = build_index::color_dbg(*this, config, n_threads);
    m_color_matrix.swap(tmp_color_matrix);

    //delete temporal files
    sdsl::util::delete_all_files(config.file_map);

    // size_type count = 1;
    // for (size_type i = 0; i < tot_nodes(); ++i) {
    //     std::cout << "id: " << i;
    //     if (m_color_matrix.is_colored(i)) {
    //         std::cout << ", col: " << count << ", label: " << node2string(i);
    //         count++;
    //     } else std::cout << ", col: " << " - " << ", label: " << node2string(i);
    //     std::cout << std::endl;
    // }

}

dbg_boss::dbg_boss(size_t K, size_t n_threads, std::string dir, std::string name) : m_n_solid_nodes(0),
                                                                                    m_k(K){

    cache_config config(true, dir, name);

    load_from_file(m_edge_bwt, cache_file_name("ebwt_wt_huff",config));
    load_from_file(m_f_array, cache_file_name("f_array",config));

    load_from_file(m_node_marks, cache_file_name("node_marks", config));
    m_node_marks_ss.set_vector(&m_node_marks);
    m_node_marks_rs.set_vector(&m_node_marks);

    load_from_file(m_solid_nodes, cache_file_name("solid_nodes", config));
    m_solid_nodes_rs.set_vector(&m_solid_nodes);
    m_solid_nodes_ss.set_vector(&m_solid_nodes);
    m_n_solid_nodes = m_solid_nodes_rs(m_solid_nodes.size());

    auto tmp_color_matrix = build_index::color_dbg(*this, config, n_threads);
    m_color_matrix.swap(tmp_color_matrix);

    //delete temporal files
    sdsl::util::delete_all_files(config.file_map);
}

dbg_boss::size_type dbg_boss::serialize(std::ostream &out, structure_tree_node *v, std::string name) const{

    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_t written_bytes = 0;
    written_bytes += m_f_array.serialize(out,child, "f_array");
    written_bytes += m_edge_bwt.serialize(out, child, "edge_bwt");
    written_bytes += m_node_marks.serialize(out, child, "node_marks");
    written_bytes += m_node_marks_rs.serialize(out, child, "node_marks_rs");
    written_bytes += m_node_marks_ss.serialize(out, child, "node_marks_ss");
    written_bytes += m_solid_nodes.serialize(out, child, "m_solid_nodes");
    written_bytes += m_solid_nodes_ss.serialize(out, child, "m_solid_nodes_ss");
    written_bytes += m_solid_nodes_rs.serialize(out, child, "m_solid_nodes_rs");
    written_bytes += m_color_matrix.serialize(out, child, "m_color_matrix");
    written_bytes += write_member(m_k, out, child, "k");
    written_bytes += write_member(m_n_solid_nodes, out, child, "m_n_solid_nodes");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void dbg_boss::load(std::istream& in){
    m_f_array.load(in);
    m_edge_bwt.load(in);
    m_node_marks.load(in);
    m_node_marks_rs.load(in, &m_node_marks);
    m_node_marks_ss.load(in, &m_node_marks);
    m_solid_nodes.load(in);
    m_solid_nodes_ss.load(in, &m_solid_nodes);
    m_solid_nodes_rs.load(in, &m_solid_nodes);
    m_color_matrix.load(in);
    read_member(m_k, in);
    read_member(m_n_solid_nodes, in);
}

dbg_boss::dbg_boss(): m_n_solid_nodes(0),
                      m_k(0){}


dbg_boss::size_type dbg_boss::rev_comp(dbg_boss::size_type v) const {

    label_t label = node_label(v);
    label_t rev_comp_label(m_k-1);

    for(size_t i=0,j=label.size()-1;i<label.size();i++,j--) rev_comp_label[i] = dna_alphabet::comp2rev[label[j]];

    //backward search return the range in the bwt not the label id!!
    std::vector<bs_t> bs_res = backward_search(rev_comp_label, 0);

    assert(!bs_res.empty());
    return bs_res[0].r_start;
}

dbg_boss::label_t dbg_boss::node_label(dbg_boss::size_type v) const {

    label_t label(m_k-1);

    for(size_t i=0;i<m_k-1;i++){
        if(v==0){
            label[i] = 1; //only dollars
            continue;
        }else {
            label[i] = m_f_array[v];
        }
        size_t rank = v - m_f_array.C[label[i]]+1;
        v = m_node_marks_rs.rank(m_edge_bwt.select(rank, label[i]<<1U));
    }
    return label;
}

std::vector<dbg_boss::bs_t> dbg_boss::backward_search(std::string query , uint8_t mismatches,
                                                      std::vector<edges_t::value_type> &syms,
                                                      std::vector<edges_t::size_type> &r_b,
                                                      std::vector<edges_t::size_type> &r_a) const {

    label_t tmp(query.size());
    for(size_t i=0;i<query.size();i++) tmp[i] = dna_alphabet::char2comp[query[i]];
    std::reverse(tmp.begin(), tmp.end());
    return backward_search(tmp, mismatches, syms, r_b, r_a);
}

std::vector<dbg_boss::bs_t> dbg_boss::backward_search(dbg_boss::label_t &query, uint8_t mismatches,
                                                      std::vector<edges_t::value_type> &syms,
                                                      std::vector<edges_t::size_type> &r_b,
                                                      std::vector<edges_t::size_type> &r_a) const {

    std::stack<bs_step_t> bs_stack;
    std::vector<bs_t>res;
    uint8_t mm, symbol;

    size_type k_start, k_end, n_symbols;

    bs_step_t root = {0, m_edge_bwt.size()-1, 0, query.size()};
    bs_stack.push(root);

    while(!bs_stack.empty()){

        bs_step_t bs_step = bs_stack.top();
        bs_stack.pop();

        if(bs_step.last_index>0) {

            m_edge_bwt.interval_symbols(bs_step.r_start, bs_step.r_end+1,
                                        n_symbols, syms, r_b, r_a);

            for(size_t i=0;i<n_symbols;i++){

                symbol =  syms[i]>>1U;
                if((syms[i] & 1U)) continue;

                mm = bs_step.mm + (symbol != query[bs_step.last_index - 1]);

                if (mm <= mismatches) {
                    k_start = m_f_array.C[symbol] + r_b[i];
                    k_end = m_f_array.C[symbol] + r_a[i];

                    bs_stack.push({m_node_marks_ss.select(k_start) + 1,
                                   m_node_marks_ss.select(k_end), mm,
                                   bs_step.last_index-1});
                }
            }
        }else {
            res.emplace_back(m_node_marks_rs(bs_step.r_start),
                             m_node_marks_rs(bs_step.r_end),
                             bs_step.mm);
        }
    }

    return res;
}

std::string dbg_boss::node2string(dbg_boss::size_type v, bool rev) const {
    label_t label = node_label(v);
    std::string label_str;
    for (unsigned char l : label) label_str.push_back(dna_alphabet::comp2char[l]);
    if(rev) std::reverse(label_str.begin(), label_str.end());
    return label_str;
}

dbg_boss::size_type dbg_boss::string2node(std::string string) const {
    assert(string.size()==(m_k-1));
    label_t label;
    std::copy(string.begin(), string.end(), std::back_inserter(label));
    std::reverse(label.begin(), label.end());
    for(unsigned char & i : label) i = dna_alphabet::char2comp[i];
    std::vector<bs_t> bs = backward_search(label,0);
    if(!bs.empty()) return bs[0].r_start;
    return 0;
}

dbg_boss::size_type dbg_boss::string2node(const char * string, size_t size) const {
    assert(size==(m_k-1));
    label_t label;
    std::copy(string, string+size, std::back_inserter(label));
    std::reverse(label.begin(), label.end());
    for(unsigned char & i : label) i = dna_alphabet::char2comp[i];
    std::vector<bs_t> bs = backward_search(label,0);
    if(!bs.empty()) return bs[0].r_start;
    return 0;
}
