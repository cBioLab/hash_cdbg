#include "build_index.hpp"
#include <bitset>
#include <sdsl/construct.hpp>
#include <ctime>
#include <cstdint>  // for me

std::mutex mtx;

inline uint64_t hash(uint64_t val1, uint64_t val2) {
    uint64_t lhs = std::hash<uint64_t>()(val1), rhs = std::hash<uint64_t>()(val2);
    uint64_t seed = lhs;
    seed ^= lhs + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);
    seed ^= rhs + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);

    return seed;
}

// uint64_t get_raw() {
//     u ^= (u << 5); u ^= (u >> 15); u ^= (u << 27);
//     w += UINT64_C(0x61C8864680B583EB);
//     return u+(w^(w>>27));
// }

inline void print_max_rss(std::string s) {
    // getrusage関数による最大消費メモリ量のチェック
    int chk;
    struct rusage usage;
    chk = getrusage(RUSAGE_SELF, &usage);
    if(chk != 0){
        printf("error\n");
        exit(-1);
    }
    // このプロセスが消費した最大メモリ領域
    printf("%sMax RSS = %lf MB\n", s.c_str(), usage.ru_maxrss / 1024.0);
    fflush(stdout);
}

void build_index::build_SA_BWT_LCP(sdsl::cache_config &config) {
    // construct SA
    sdsl::construct_sa<8>(config);
    register_cache_file(sdsl::conf::KEY_SA, config);

    // construct BWT
    sdsl::construct_bwt<8>(config);
    register_cache_file(sdsl::conf::KEY_BWT, config);

    // construct LCP
    sdsl::construct_lcp_semi_extern_PHI(config);
    register_cache_file(sdsl::conf::KEY_LCP, config);
}

void build_index::build_edge_BWT(sdsl::cache_config &config, size_t K) {

    bool flag_symbols;
    std::bitset<6> kmer_symbols;
    std::bitset<6> flagged_symbols=false;
    sdsl::bit_vector node_marks;
    sdsl::int_vector_buffer<8> bwt(cache_file_name(sdsl::conf::KEY_BWT, config));
    sdsl::int_vector_buffer<8> ebwt(cache_file_name("ebwt", config), std::ios::out);
    sdsl::int_vector_buffer<> klcp(cache_file_name("klcp", config));
    size_t bwt_pos=0;

    //this is for marking dollar-prefixed kmers
    sdsl::sd_vector<> dollar_bv;
    sdsl::bit_vector solid_marks;
    load_from_file(dollar_bv, cache_file_name("sp_dollar_bv", config));
    sdsl::sd_vector<>::rank_1_type dr(&dollar_bv);
    sdsl::int_vector_buffer<> sa(cache_file_name(sdsl::conf::KEY_SA, config));
    size_t node_pos=0;

    //
    size_t n_dollars = dr(dollar_bv.size());

    //TODO just testing
    //int_vector_buffer<8> text(cache_file_name(conf::KEY_TEXT, config));
    //

    sdsl::util::assign(node_marks, sdsl::int_vector<1>(bwt.size(), 1));
    sdsl::util::assign(solid_marks, sdsl::int_vector<1>(bwt.size(), 0));

    //TODO just testing
    /*for(size_t i=0;i<20;i++){
        std::cout<<klcp[i]<<" ";
        for(size_t j=sa[i];j<std::min(sa[i]+50, text.size());j++){
            std::cout<<dna_alphabet::comp2char[text[j]];
        }
        std::cout<<""<<std::endl;
    }*/
    //

    //set the dummy node
    kmer_symbols[1]=true;
    for(size_t i=0;i<n_dollars;i++){
        if(!kmer_symbols[bwt[i]]) kmer_symbols[bwt[i]]=true;
        //62=111110
        if(kmer_symbols==62)break;
    }

    for(size_t j=1;j<kmer_symbols.size();j++) {
        if (kmer_symbols[j]) {
            ebwt[bwt_pos] = (j << 1U);
            bwt_pos++;
        }
    }
    node_marks[bwt_pos-1] = false;
    node_pos++;
    kmer_symbols.reset();
    //

    for(size_t i=2;i<=klcp.size();i++){

        //TODO just checking
        //bool first=true;

        if((i==klcp.size() || klcp[i]<K) && kmer_symbols!=0){

            //this is for marking solid nodes
            if(dr.rank(std::min(sa[i-1]+K, sa.size()))-dr.rank(sa[i-1])==0){
                solid_marks[node_pos]=true;
            }


            if(i<=(n_dollars+1)){

                ebwt[bwt_pos] = 3;
                bwt_pos++;

                //TODO just testing
                /*if(color_marks[node_pos]) {
                    std::cout << epos << "\t" << (ebwt[bwt_pos - 1] & 1UL)
                              << (dna_alphabet::comp2char[ebwt[bwt_pos - 1] >> 1U]) << "\t"
                              << i << "\t";

                    epos++;
                    bool dollar_start = false;
                    for (size_t i1 = 0; i1 < K; i1++) {
                        if (((sa[i - 1] + i1) >= text.size() || (text[sa[i - 1] + i1]) == 1) && text[sa[i - 1]] != 1) {
                            dollar_start = true;
                        }

                        if (!dollar_start) {
                            std::cout << dna_alphabet::comp2char[text[sa[i - 1] + i1]];
                        } else {
                            std::cout << "$";
                        }
                    }
                    std::cout << "\t" << solid_marks[node_pos] << "\t" << color_marks[node_pos]<<"\t"<<node_pos;
                    kpos++;
                    std::cout << "" << std::endl;
                }*/
                //

            }else {
                //store the symbols of the previous kmer
                for(size_t j=1;j<kmer_symbols.size();j++){

                    if(kmer_symbols[j]){
                        ebwt[bwt_pos] = (j << 1U) | flagged_symbols[j];
                        bwt_pos++;

                        //TODO just checking
                        /*if((node_pos>=50810 && node_pos<50820) || node_pos==18788) {
                            std::cout << bwt_pos-1 << "\t" << (ebwt[bwt_pos - 1] & 1UL)
                                      << (dna_alphabet::comp2char[ebwt[bwt_pos - 1] >> 1U]) << "\t"
                                      << i << "\t";
                            epos++;
                            if (first) {
                                bool dollar_start = false;
                                for (size_t i1 = 0; i1 < K; i1++) {

                                    if (((sa[i - 1] + i1) >= text.size() || (text[sa[i - 1] + i1]) == 1) &&
                                        text[sa[i - 1]] != 1) {
                                        dollar_start = true;
                                    }

                                    if (!dollar_start) {
                                        std::cout << dna_alphabet::comp2char[text[sa[i - 1] + i1]];
                                    } else {
                                        std::cout << "$";
                                    }
                                }
                                std::cout << "\t" << solid_marks[node_pos]<<"\t"<<color_marks[node_pos]<<"\t"<<node_pos;
                                kpos++;
                            }
                            first = false;
                            std::cout << "" << std::endl;
                        }*/
                        //
                    }
                }
            }

            node_marks[bwt_pos-1] = false;
            node_pos++;
            if(i==klcp.size()) break;

            //setup for the next Kmer;
            flag_symbols = klcp[i] == (K - 1);
            if (!flag_symbols) {
                flagged_symbols.reset();
            } else {
                flagged_symbols |= kmer_symbols;
            }
            kmer_symbols.reset();
        }
        kmer_symbols.set(bwt[i]);
    }

    node_marks.resize(bwt_pos);
    solid_marks.resize(node_pos);

    ebwt.close();
    klcp.close();

    sdsl::wt_huff<> tmp_wt;
    sdsl::construct(tmp_wt, cache_file_name("ebwt", config));

    sdsl::int_vector<64> tmp_C;
    tmp_C.resize(dna_alphabet::sigma+1);

    //compute the rank of every symbol
    tmp_C[0]=0;// # symbol
    for(uint8_t i=1;i<dna_alphabet::sigma;i++) tmp_C[i] = tmp_wt.rank(tmp_wt.size(),(i<<1U));

    uint64_t tmp1, tmp2;
    tmp2 = 0;
    for(uint8_t i=1;i<=dna_alphabet::sigma;i++){
        tmp1 = tmp_C[i];
        tmp_C[i] = tmp_C[i-1] + tmp2;
        tmp2 = tmp1;
    }

    /*size_t tmp_cont=0;
    for(auto && color_mark : color_marks){
        if(color_mark){
            tmp_cont++;
        }
    }
    std::cout<<"estimated what?:"<<tmp_cont<<std::endl;*/

    //wavelet three with the edges
    store_to_file(tmp_wt, cache_file_name("ebwt_wt_huff", config));

    //first column (F array in the FM-index)
    f_array_t tmp_alph(tmp_C);
    store_to_file(tmp_alph, cache_file_name("f_array", config));

    //bit vector marking solid nodes
    sdsl::rrr_vector<63> rrr_solid_marks(solid_marks);
    store_to_file(rrr_solid_marks, cache_file_name("solid_nodes", config));

    //bit vector marking edges
    sdsl::rrr_vector<63> tmp_node_marks(node_marks);
    store_to_file(tmp_node_marks, cache_file_name("node_marks", config));

    //remove unnecessary files
    remove(cache_file_name(sdsl::conf::KEY_SA, config));
    remove(cache_file_name(sdsl::conf::KEY_BWT, config));
    remove(cache_file_name("ebwt", config));
    remove(cache_file_name("klcp", config));
    remove(cache_file_name("sp_dollar_bv", config));

    register_cache_file("solid_nodes", config);
    register_cache_file("node_marks", config);
    register_cache_file("ebwt_wt_huff", config);
    register_cache_file("f_array", config);
}

void build_index::build_KLCP(sdsl::cache_config &config, size_t K) {

    size_t n_dollars;
    sdsl::int_vector_buffer<> lcp(cache_file_name(sdsl::conf::KEY_LCP, config));
    sdsl::int_vector_buffer<> klcp(cache_file_name("klcp", config), std::ios::out, 1000000);

    sdsl::int_vector_buffer<> sa(cache_file_name(sdsl::conf::KEY_SA, config));
    sdsl::sd_vector<> dollar_bv;
    load_from_file(dollar_bv, cache_file_name("sp_dollar_bv", config));
    sdsl::sd_vector<>::rank_1_type dr(&dollar_bv);
    sdsl::sd_vector<>::select_1_type ds(&dollar_bv);
    n_dollars = dr.rank(dollar_bv.size());

    klcp[0] = 0;

    //dummy node
    //TODO K cannot be greater than the smallest read!
    for(size_t i=1;i<=n_dollars;i++){
        if(lcp[i]>K){
            klcp[i] = K;
        }else{
            klcp[i] = lcp[i];
        }
    }

    for(size_t i=(n_dollars+1);i<lcp.size();i++){
        //modify the LCP
        if(lcp[i]>K || lcp[i]>=(ds.select(dr.rank(sa[i]+1)+1)-sa[i]+1)){
            klcp[i]=K;
        }else{
            klcp[i] = lcp[i];
        }
    }

    lcp.close();
    klcp.close();

    remove(cache_file_name(sdsl::conf::KEY_LCP, config));
}

void build_index::build_boss(std::string input_file, sdsl::cache_config &config, size_t K) {
    fastx_parser::preproc_reads(input_file, config, 1);
    build_SA_BWT_LCP(config);
    build_KLCP(config, K - 1);
    build_edge_BWT(config, K - 1);
}

build_index::color_matrix_t build_index::color_dbg(dbg_boss &dbg_index,
                                                   cache_config &config,
                                                   size_t n_threads) {
    std::cout << "test hash function" << std::endl;

    size_t elm_per_thread, rem, start, end;

    matrix_skeleton m_skeleton;

    print_max_rss("[color_dbg] ");

    std::cout<<"Computing the skeleton of the color matrix"<<std::endl;
    build_matrix_skeleton(dbg_index, m_skeleton, n_threads, config);

    print_max_rss("[color_dbg] ");

    //coloring the nodes
    std::cout<<"Filling the matrix"<<std::endl;
    std::vector<std::thread> threads;
    text_t text(cache_file_name(conf::KEY_TEXT, config));
    elm_per_thread = text.size()/n_threads;
    rem = text.size() % n_threads;
    text.close();

    for(size_t i=0;i<n_threads;i++){
        start = i*elm_per_thread;
        end = (i+1)*elm_per_thread;
        if(i==n_threads-1) end +=rem-1;
        threads.emplace_back(std::thread(color_dbg_internal,
                                         std::ref(dbg_index),
                                         start,end,
                                         std::ref(m_skeleton),
                                         std::ref(config)));
    }
    for (auto & it : threads) it.join();
    //

    print_max_rss("[color_dbg] ");

    //store the compressed version of the color matrix
    std::cout<<"Compacting the matrix"<<std::endl;
/*
    size_type n_rows = m_skeleton.colored_rows_rs(m_skeleton.colored_rows.size());

    for(size_type j=1;j<=n_rows;j++){

        start = m_skeleton.colored_cells_ss(j);
        end = m_skeleton.colored_cells_ss(j+1)-1;
        std::cout << "node_j: " << j << ", color start: " << start << ", end: " << end << std::endl;

        // duple color チェックどうする？

        // std::vector<uint64_t> node_colors(end-start+1);  // 16 to 64
        // std::cout << "node j " << j << ", colors:";  // for debug
        for(size_t i=start,k=0;i<=end;i++,k++) {
            // std::cout << " " << m_skeleton.cell_values[i];  // for debug
            // for me
            // if (m_skeleton.cell_values[i] == 0) {
            //     k--;
            //     continue;
            // }

            std::cout << "  " << m_skeleton.cell_values[i] << std::endl;
            // if (m_skeleton.cell_values[i] == 0) {
                // k--;
                // std::cout << "color 0!" << std::endl;
            // } else {
            //     node_colors[k] = m_skeleton.cell_values[i];
            // }
            // assert(node_colors[k]!=0);
        }
        // std::cout << std::endl;  // for debug

        // std::sort(node_colors.begin(), node_colors.end());

        // for(size_t i=start,k=0;i<=end;i++,k++){
        //     if (node_colors[k] == 0) {
        //         std::cout << "nodecolor 0!" << std::endl;
        //         i--;
        //         continue;
        //     }
        //     if(i==start){
        //         m_skeleton.cell_values[i] = node_colors[k];
        //     }else{
        //         m_skeleton.cell_values[i] = node_colors[k]-node_colors[k-1];
        //         // for me
        //         // assert((node_colors[k]-node_colors[k-1])!=0);
        //         if ((node_colors[k] - node_colors[k - 1] == 0)) {
        //             std::cout << j << ": " << node_colors[k] << ", " << node_colors[k - 1] << " duple color!" << std::endl;
        //         }
        //     }
        // }
    }
*/
    // color_matrix_t cm(m_skeleton, n_threads);
    color_matrix_t cm(m_skeleton, std::ref(dbg_index), n_threads);
    return cm;
}

void build_index::color_dbg_internal(dbg_boss &dbg_index,
                                     long long int start,
                                     long long int end,
                                     matrix_skeleton& m_skeleton,
                                     sdsl::cache_config &config) {

    text_t text(cache_file_name(conf::KEY_TEXT, config));
    size_type tmp_node, node_row, row_pos, last_row_pos, tmp1, tmp2, middle;

    dbg_boss::label_t label(dbg_index.k-1);


    std::vector<uint8_t> syms(dna_alphabet::sigma*2);
    std::vector<size_type> r_b(dna_alphabet::sigma*2);
    std::vector<size_type> r_a(dna_alphabet::sigma*2);

    //I am using a buffer which is big enough
    std::vector<std::pair<size_type, size_type> > col_ranges(10000);

    uint64_t new_color;  // 32 to 64
    long long int text_pos;

    //PARALLEL
    //adjust start and end
    while(text[start]!=1) start--;
    while(text[end]!=1) end--;

    // for me
    std::map<uint64_t, bool> used_colors_bv;
    std::vector<size_type> order;
    // std::vector<size_type> nodes;  // for debug
    // for debug
    // long long int ambdiv = (end - start) / 10;

    text_pos=end;
    while(text_pos>start){

        for(size_t j=dbg_index.k-1;j-->0;){
            label[j] = text[text_pos];
            text_pos--;
        }

        //the node labelled with the start of the read
        tmp_node = dbg_index.backward_search(label,
                                             0, syms,
                                             r_b, r_a)[0].r_start;

        assert(tmp_node!=0 && m_skeleton.colored_rows[tmp_node]);

        //the starting dBG node of the read
        order.emplace_back(tmp_node);  // for me
        // nodes.emplace_back(tmp_node);  // for debug

        // std::cout << "\n== check new seq ==" << std::endl;  // for debug
        // std::cout << "check node " << tmp_node << ", " << dbg_index.node2string(tmp_node) << std::endl;  // for debug

        while(text[text_pos+1]!=1){
            tmp_node = dbg_index.outgoing(tmp_node,
                                          text[text_pos],
                                          true, syms, r_b, r_a);

            // std::cout << "check node " << tmp_node << ", " << dbg_index.node2string(tmp_node) << std::endl;  // for debug

            if(m_skeleton.colored_rows[tmp_node]){
                order.emplace_back(tmp_node);
                // nodes.emplace_back(tmp_node);  // for debug
            }

            text_pos--;
        }

        size_t cont=0;
        for (auto const &node : order) {
            node_row = m_skeleton.colored_rows_rs(node);
            col_ranges[cont] = std::make_pair(m_skeleton.colored_cells_ss(node_row + 1),
                                              m_skeleton.colored_cells_ss(node_row + 2) - 1);
            cont++;
        }
        //END PARALLEL


        //SEQUENTIAL
        try {
            //retrieve the get_colors that we cannot use for the current read
            std::lock_guard<std::mutex> lck(mtx);

            //find all used colors
            row_pos = col_ranges[0].first;
            last_row_pos = col_ranges[0].second;
            uint64_t color;
            while(row_pos <= last_row_pos){
                color = m_skeleton.cell_values[row_pos];
                if(color==0)break;
                used_colors_bv[color] = true;
                row_pos++;
            }

            //find the next available color
            // new_color = hash(order[0], (uint64_t)text_pos);
            new_color = 1;
            // std::cout << std::endl;
            // std::cout << "node: " << nodes[0] << std::endl;
            while (used_colors_bv[new_color]) {
                // new_color = hash(order[0], new_color);
                new_color++;
                if (new_color == 0) {
                    new_color = 1;
                }
            }
            // std::cout << ", new color: " << new_color << std::endl;
            assert(new_color > 0 && new_color < UINT64_MAX - 1);

            //assign the new color to the nodes that need to be colored
            // for me
            size_t color_idx = 0;
            for (auto const &node : order) {
                // std::cout << "paint node " << nodes[color_idx] << ", col: " << new_color << ", " << dbg_index.node2string(nodes[color_idx]) << std::endl;  // for debug
                row_pos = col_ranges[color_idx].first;
                last_row_pos = col_ranges[color_idx].second;
                color_idx++;
                // std::cout << "num: " << last_row_pos - row_pos << std::endl;  // for debug

                //binary search the next available position
                tmp1 = row_pos;
                tmp2 = last_row_pos;
                while (tmp1 <= tmp2) {
                    middle = tmp1 + (tmp2 - tmp1) / 2;

                    // Check if x is present at mid
                    if (m_skeleton.cell_values[middle] == 0 &&
                         (middle==tmp1 || m_skeleton.cell_values[middle-1]!=0 )){
                        row_pos=middle;
                        break;
                    }

                    // If x greater, ignore left half
                    if (m_skeleton.cell_values[middle] != 0){
                        tmp1 = middle + 1;
                    }else{
                        tmp2 = middle - 1;
                    }
                }
                assert(row_pos<=last_row_pos);
                m_skeleton.cell_values[row_pos]=new_color;

                // std::cout << "change color, node: " << nodes[color_idx] << ", hash: " << node << " and " << new_color;  // for debug
                new_color = hash(node, new_color);
                // std::cout << ", new: " << new_color << std::endl;  // for debug
            }
            //END SEQUENTIAL
        }catch (std::logic_error&){
            std::cout<<"some error in the threads"<<std::endl;
        }

        text_pos++;

        // for me
        used_colors_bv.clear();
        order.clear();
        // nodes.clear();  // for debug

        // for debug
        // if ((text_pos - start) % ambdiv == 0) {
        //     std::cout << "start: " << start << " of " << (text_pos - start) / ambdiv << std::endl;
        // }
    }
}

void build_index::estimate_color_rows(dbg_boss &dbg_index,
                                      size_type start,
                                      size_type end,
                                      sdsl::bit_vector &color_marks) {

    size_t ind;

    for(size_type i=start; i<=end;i++){

        if(i==0) continue;

        if(!dbg_index.solid_nodes[i]){
            if(dbg_index.solid_nodes[dbg_index.incomming(i,1)]){
                color_marks[i] = true;
            }else if(dbg_index.solid_nodes[dbg_index.outgoing(i,1)]){
                color_marks[i] = true;
            }
        }else{
            ind = dbg_index.indegree(i);
            for(size_t j=1;j<=ind;j++){
                if(dbg_index.outdegree(dbg_index.incomming(i,j))>1){
                    color_marks[i] = true;
                }
            }
        }
    }
}

void build_index::estimate_color_cells(dbg_boss &dbg_index,
                                       long long int start,
                                       long long int end,
                                       sdsl::cache_config &config,
                                       sdsl::rrr_vector<63> &color_marks,
                                       sdsl::rrr_vector<63>::rank_1_type &color_marks_rs,
                                    //    std::map<build_index::size_type, std::map<build_index::size_type, bool> > &n_nodes,
                                       std::map<build_index::size_type, uint64_t> &c_nodes) {  // 32 to 64

    text_t text(cache_file_name(conf::KEY_TEXT, config));

    //adjust start and end
    while(text[start]!=1) start--;
    while(text[end]!=1) end--;
    long long int text_pos;

    //temporal files
    dbg_boss::label_t label(dbg_index.k-1);
    std::vector<uint8_t> syms(dna_alphabet::sigma*2);
    std::vector<size_type> r_b(dna_alphabet::sigma*2);
    std::vector<size_type> r_a(dna_alphabet::sigma*2);
    size_type tmp_node;

    std::map<size_type, long long> kmers_in_read;  // bool to long long
    text_pos=end;

    // std::map<size_type, std::map<size_type, bool> > check_node;
    // size_t ind, outd;
    // size_type prev_node, nei_node;

    while(text_pos>start) {

        for(size_t j=dbg_index.k-1;j-->0;){
            label[j] = text[text_pos];
            text_pos--;
        }
        tmp_node = dbg_index.backward_search(label,
                                             0, syms,
                                             r_b, r_a)[0].r_start;

        kmers_in_read[color_marks_rs(tmp_node)]++;  // true to ++

        while(text[text_pos+1]!=1) {

            tmp_node = dbg_index.outgoing(tmp_node,
                                          text[text_pos],
                                          true, syms, r_b, r_a);
            if(color_marks[tmp_node]){
                //for me
                // kmers_in_read[color_marks_rs(tmp_node)]=true;
                kmers_in_read[color_marks_rs(tmp_node)]++;

                // ind = dbg_index.indegree(tmp_node);
                // for (size_t i = 1; i <= ind; ++i) {
                //     prev_node = dbg_index.incomming(tmp_node, i, false);
                //     outd = dbg_index.outdegree(prev_node);
                //     if(outd > 1) {
                //         for (size_t j = 1; j <= outd; ++j) {
                //             nei_node = dbg_index.outgoing(prev_node, j, false);
                //             if (nei_node == tmp_node) continue;
                //             // check_node[color_marks_rs(tmp_node)][color_marks_rs(nei_node)] = true;
                //             n_nodes[color_marks_rs(tmp_node)][color_marks_rs(nei_node)] = true;
                //         }
                //     }
                // }
            }

            text_pos--;
        }

        for(auto const &kmer : kmers_in_read){
            // for me
            // c_nodes[kmer.first]++;
            c_nodes[kmer.first] += kmer.second;
        }

        // for (auto const &node : check_node) {
        //     for (auto const &nei : node.second) {
        //         n_nodes[node.first][nei.first] = true;
        //     }
        // }

        kmers_in_read.clear();
        // check_node.clear();
        text_pos++;
    }
}

void build_index::build_matrix_skeleton(dbg_boss &dbg_index,
                                        matrix_skeleton &skeleton,
                                        size_t n_threads,
                                        sdsl::cache_config &config) {


    size_t elm_per_thread, rem, start, end;

    //marking the nodes that need to be colored can be done in parallel
    elm_per_thread = dbg_index.tot_nodes()/n_threads;
    rem = dbg_index.tot_nodes() % n_threads;

    print_max_rss("[bui_m_ske] ");

    {
        bit_vector colored_rows_bv;
        sdsl::util::assign(colored_rows_bv, sdsl::int_vector<1>(dbg_index.tot_nodes(), 0));
        std::vector<std::thread> col_marks_threads;
        for (size_t i = 0; i < n_threads; i++) {
            start = i * elm_per_thread;
            end = (i + 1) * elm_per_thread - 1;
            if (i == n_threads - 1) end += rem;
            col_marks_threads.emplace_back(std::thread(estimate_color_rows,
                                                       std::ref(dbg_index),
                                                       start,
                                                       end,
                                                       std::ref(colored_rows_bv)));
        }
        for (auto &it : col_marks_threads) it.join();
        sdsl::rrr_vector<63> colored_rows(colored_rows_bv);
        skeleton.colored_rows.swap(colored_rows);
    }
    skeleton.colored_rows_rs.set_vector(&skeleton.colored_rows);

    print_max_rss("[bui_m_ske] ");

    {
        text_t text(cache_file_name(conf::KEY_TEXT, config));
        elm_per_thread = text.size() / n_threads;
        rem = text.size() % n_threads;
        text.close();

        std::map<size_type, uint64_t> all_color_map;  // 32 to 64
        size_type tot_color_cells = 0;
        {
            //building the array that will contain the uncompressed colors
            std::vector<std::thread> col_est_threads;
            std::vector<std::map<size_type, uint64_t>> color_maps(n_threads);  // 32 to 64
            // std::vector<std::map<size_type, std::map<size_type, bool> > > neighbor_maps(n_threads);
            for (size_t i = 0; i < n_threads; i++) {
                start = i * elm_per_thread;
                end = (i + 1) * elm_per_thread;
                if (i == n_threads - 1) end += rem - 1;
                col_est_threads.emplace_back(std::thread(estimate_color_cells,
                                                         std::ref(dbg_index),
                                                         start, end,
                                                         std::ref(config),
                                                         std::ref(skeleton.colored_rows),
                                                         std::ref(skeleton.colored_rows_rs),
                                                        //  std::ref(neighbor_maps[i]),
                                                         std::ref(color_maps[i])));
            }
            for (auto &it : col_est_threads) it.join();

            print_max_rss("[bui_m_ske] ");

            //count the number of colors
            for (size_t i = 0; i < n_threads; i++) {
                for (auto const &node: color_maps[i]) {
                    all_color_map[node.first] += node.second;
                    tot_color_cells += node.second;
                }
                color_maps[i].clear();

                // for (auto const &node : neighbor_maps[i]) {
                //     for (auto const &nei : node.second) {
                //         skeleton.neighbor_nodes[node.first][nei.first] = true;
                //     }
                // }
            }
            color_maps.clear();

            //TODO testing
            /*std::vector<std::thread> amb_threads;
            size_t n_amb=0;
            for (size_t i = 0; i < n_threads; i++) {
                start = i * elm_per_thread;
                end = (i + 1) * elm_per_thread;
                if (i == n_threads - 1) end += rem - 1;
                amb_threads.emplace_back(std::thread(estimate_amb_seqs,
                                                     std::ref(dbg_index),
                                                     start, end,
                                                     std::ref(config),
                                                     std::ref(skeleton.colored_rows),
                                                     std::ref(skeleton.colored_rows_rs),
                                                     std::ref(n_amb)));
            }
            for (auto &it : amb_threads) it.join();
            std::cout<<"there are "<<n_amb<<" ambiguous reads"<<std::endl;*/
            //
        }

        print_max_rss("[bui_m_ske] ");

        bit_vector colored_cells;
        sdsl::util::assign(colored_cells, sdsl::int_vector<1>(tot_color_cells + 1, 0));

        size_type pos = 0;
        colored_cells[pos] = true;

        for (auto const &node: all_color_map) {
            pos += node.second;
            colored_cells[pos] = true;
        }
        all_color_map.clear();

        sdsl::rrr_vector<63> comp_colored_cells(colored_cells);
        skeleton.colored_cells.swap(comp_colored_cells);
        skeleton.colored_cells_ss.set_vector(&skeleton.colored_cells);
        //
    }

    print_max_rss("[bui_m_ske] ");

    // sdsl::util::assign(skeleton.cell_values, sdsl::int_vector<64>(skeleton.colored_cells.size()-1, 0));  // 16 to 64
    // skeleton.cell_values.resize(skeleton.colored_cells.size()-1);
    // for (auto &val : skeleton.cell_values) {
    //     if (val != 0) val = 0;
    // }
    sdsl::int_vector<64> tmp_cell_values(skeleton.colored_cells.size()-1, 0);
    skeleton.cell_values.swap(tmp_cell_values);

    print_max_rss("[bui_m_ske] ");
}

void
build_index::estimate_amb_seqs(dbg_boss &dbg_index, long long int start, long long int end, sdsl::cache_config &config,
                               sdsl::rrr_vector<63> &color_marks, sdsl::rrr_vector<63>::rank_1_type &color_marks_rs,
                               size_t& tot_amb) {

    text_t text(cache_file_name(conf::KEY_TEXT, config));

    //adjust start and end
    while(text[start]!=1) start--;
    while(text[end]!=1) end--;
    long long int text_pos;
    bool is_amb;
    size_t n_amb=0, outd;

    //temporal files
    dbg_boss::label_t label(dbg_index.k-1);
    std::vector<uint8_t> syms(dna_alphabet::sigma*2);
    std::vector<size_type> r_b(dna_alphabet::sigma*2);
    std::vector<size_type> r_a(dna_alphabet::sigma*2);
    size_type tmp_node, tmp_node2, neighbor_node;

    std::map<size_type, bool> kmers_in_read;
    std::map<size_type, bool> neighbor_colored;
    text_pos=end;

    while(text_pos>start) {

        for(size_t j=dbg_index.k-1;j-->0;){
            label[j] = text[text_pos];
            text_pos--;
        }

        tmp_node = dbg_index.backward_search(label,
                                             0, syms,
                                             r_b, r_a)[0].r_start;
        outd = dbg_index.outdegree(tmp_node);

        if(outd>1){
            tmp_node2 = dbg_index.outgoing(tmp_node,
                                           text[text_pos],
                                           true, syms, r_b, r_a);
            for(size_t i=1;i<=outd;i++){
                neighbor_node = dbg_index.outgoing(tmp_node, i, false);
                if(neighbor_node!=tmp_node2 && color_marks[neighbor_node]){
                    neighbor_colored[color_marks_rs(neighbor_node)]=true;
                }
            }
        }

        kmers_in_read[color_marks_rs(tmp_node)]=true;
        is_amb=false;

        while(text[text_pos+1]!=1) {

            tmp_node = dbg_index.outgoing(tmp_node,
                                          text[text_pos],
                                          true, syms, r_b, r_a);
            if(color_marks[tmp_node]){

                if(kmers_in_read[color_marks_rs(tmp_node)]){
                    is_amb=true;
                }else{
                    kmers_in_read[color_marks_rs(tmp_node)]=true;
                }
            }

            outd = dbg_index.outdegree(tmp_node);
            if(outd>1 && dbg_index.solid_nodes[tmp_node]){
                tmp_node2 = dbg_index.outgoing(tmp_node,
                                              text[text_pos-1],
                                              true, syms, r_b, r_a);
                for(size_t i=1;i<=outd;i++){
                    neighbor_node = dbg_index.outgoing(tmp_node, i, false);
                    if(neighbor_node!=tmp_node2 && color_marks[neighbor_node]){
                        neighbor_colored[color_marks_rs(neighbor_node)]=true;
                    }
                }
            }
            text_pos--;
        }

        for(auto const& kmers: kmers_in_read){
            if(neighbor_colored[kmers.first]){
                is_amb = true;
                break;
            }
        }

        if(is_amb){
            n_amb++;
        }

        neighbor_colored.clear();
        kmers_in_read.clear();
        text_pos++;
    }

    try {
        //retrieve the get_colors that we cannot use for the current read
        std::lock_guard<std::mutex> lck(mtx);
        tot_amb+=n_amb;
    }catch (std::logic_error&){
        std::cout<<"some error in the threads"<<std::endl;
    }

}
