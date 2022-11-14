#include <hash_dbg/dna_alphabet.hpp>
#include <hash_dbg/boss.hpp>
#include <fstream>
#include <thread>
#include <mutex>
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>

std::mutex mtx_seqs;

inline uint64_t hash(uint64_t val1, uint64_t val2) {
    uint64_t lhs = std::hash<uint64_t>()(val1), rhs = std::hash<uint64_t>()(val2);
    uint64_t seed = lhs;
    seed ^= lhs + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);
    seed ^= rhs + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);

    return seed;
}

void compute_starters(size_t start,
                      size_t end,
                      dbg_boss& my_dbg,
                      std::unordered_map<dbg_boss::size_type, bool>& starter_hash){

    dbg_boss::size_type prev_node;

    for(size_t i=start;i<=end;i++) {
        prev_node = my_dbg.incomming(my_dbg.solid_nodes_ss(i), '$', true);
        if(prev_node!=0){
            try {
                std::lock_guard<std::mutex> lck(mtx_seqs);
                if (starter_hash.count(prev_node) == 0) {
                    starter_hash[prev_node] = true;
                }
            }catch (std::logic_error&){
                std::cout<<"some error in the threads"<<std::endl;
            }
        }
    }
}

void build_sequences(size_t start,
                     size_t end,
                     size_t& seqs_rebuilt,
                     size_t& amb_seqs,
                     size_t& fm_seqs,
                     size_t& nf_seqs,
                     std::vector<dbg_boss::size_type>& starter_nodes,
                     dbg_boss& my_dbg,
                     sdsl::csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>& fm_index,
                     std::ofstream& output_file){

    dbg_boss::size_type first_node, tmp_node, neighbor, cand_node;
    size_t outd, n_colors, n_ambigous=0, n_fm=0, n_nf=0, seq_counter=0;
    bool ambigous;
    uint8_t symbol;
    std::stringstream ss;
    uint64_t seq_color;


    size_t tmp=0;

    // for debug
    size_t ambdiv = (end - start) / 10;

    for(size_t i=start;i<=end;i++){

        first_node = starter_nodes[i];

        if(first_node!=0){  // tmp_node is the start of one or more sequences
            for (uint64_t node_color = 1; node_color <= my_dbg.get_color_num(first_node); ++node_color) {
                tmp_node = first_node;
                seq_color = node_color;
                bool flag = true;
                ambigous = false;
                std::string tmp_string, tmp_rev_string;  // strout
                tmp_string = my_dbg.node2string(tmp_node);  // strout
                for(auto & sym : tmp_string) sym = dna_alphabet::char2comp[sym];  // strout

                while (flag) {
                    if (my_dbg.is_colored(tmp_node)) {
                        seq_color = hash(tmp_node, seq_color);
                    }

                    outd = my_dbg.outdegree(tmp_node);
                    if (outd == 1) {
                        symbol = my_dbg.edge_bwt[my_dbg.get_edges(tmp_node).first] >> 1U;  // strout
                        tmp_string.push_back(symbol);  // strout
                        tmp_node = my_dbg.outgoing(tmp_node, 1);
                    } else {
                        n_colors = 0;

                        for (size_t j = 1; j <= outd; j++) {
                            neighbor = my_dbg.outgoing(tmp_node, j);
                            if (my_dbg.is_contained(neighbor, seq_color)) {
                                cand_node = neighbor;
                                symbol = my_dbg.edge_bwt[my_dbg.get_edges(tmp_node).first + j - 1] >> 1U;  // strout
                                n_colors++;
                            }
                        }

                        if (n_colors == 1) {
                            tmp_node = cand_node;
                            tmp_string.push_back(symbol);  // strout
                        } else if (n_colors == 0) {
                            ambigous = true;
                            flag = false;
                            n_nf++;
                            break;
                        } else {
                            ambigous = true;
                            flag = false;
                            n_ambigous++;
                            break;
                        }
                    }

                    if (!my_dbg.solid_nodes[tmp_node]) {
                        flag = false;
                    }
                }

                if (!ambigous) {

                    tmp_string = tmp_string.substr(1, tmp_string.size()-2);  // strout
                    seq_counter++;

                    size_t fm_count = sdsl::count(fm_index,tmp_string);  // strout

// /*  // strout
                    if(fm_count==0){
                        tmp_rev_string = tmp_string;
                        std::reverse(tmp_rev_string.begin(), tmp_rev_string.end());
                        for(auto& sym: tmp_rev_string){
                            sym = dna_alphabet::comp2rev[sym];
                        }
                        fm_count = sdsl::count(fm_index,tmp_rev_string);
                    }

                    if(fm_count==0){
                        for(auto&sym : tmp_string) sym = dna_alphabet::comp2char[sym];
                        std::cout<<tmp_string<<std::endl;
                    }
                    if (fm_count == 0) {
                        n_fm++;
                        seq_counter--;
                        std::cout << "fm_count 0!" << std::endl;
                        continue;
                    }
// */
                    for (auto& sym: tmp_string) sym = dna_alphabet::comp2char[sym];
                    ss << ">first_node_id:" << first_node << '\n';
                    ss << tmp_string << '\n';
                }
            }
        }

        // for debug
        if (ambdiv > 0 && (i - start) % ambdiv == 0) {
            std::cout << "start: " << start << " of " << (i - start) / ambdiv << ", reb: " << seq_counter << ", amb: " << n_ambigous << ", fm: " << n_fm << ", nf: " << n_nf << std::endl;
        }
    }

    size_t printed_syms = 1;
    if(printed_syms!=0){
        try {
            std::lock_guard<std::mutex> lck(mtx_seqs);
            output_file << ss.str();
        }catch (std::logic_error&){
            std::cout<<"some error in the threads"<<std::endl;
        }
        ss.clear();
    }

    try {
        std::lock_guard<std::mutex> lck(mtx_seqs);
        seqs_rebuilt+=seq_counter;
        amb_seqs += n_ambigous;
        fm_seqs += n_fm;
        nf_seqs += n_nf;
    }catch (std::logic_error&){
        std::cout<<"some error in the threads"<<std::endl;
    }
}

int main(int argc, char* argv[]) {

    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " INPUT_BOSS_INDEX INPUT_FM_INDEX N_TREADS OUTPUT_FILE" << std::endl;
        return 1;
    }

    dbg_boss dbg_index;
    sdsl::load_from_file(dbg_index, argv[1]);

    sdsl::csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;
    sdsl::load_from_file(fm_index, argv[2]);

    size_t n_threads = size_t(std::stoi(argv[3]));
    std::ofstream output_file;

    output_file.open(std::string(argv[4]) + ".fasta");

    std::vector<dbg_boss::size_type> starter_nodes;

    dbg_boss::size_type start, end;
    size_t rem = dbg_index.n_solid_nodes % n_threads;

    std::cout<<"computing the start of every sequence"<<std::endl;
    {
        dbg_boss::size_type nodes_per_thread = dbg_index.n_solid_nodes / n_threads;
        std::vector<std::thread> threads;
        std::unordered_map<dbg_boss::size_type, bool> starters_hash;
        for (size_t i = 0; i < n_threads; i++) {
            start = i * nodes_per_thread + 1;
            end = ((i + 1) * nodes_per_thread);
            if (i == n_threads - 1) end += rem;

            threads.emplace_back(std::thread(compute_starters,
                                             start,
                                             end,
                                             std::ref(dbg_index),
                                             std::ref(starters_hash)));

        }
        for (auto &it : threads) it.join();

        starter_nodes.reserve(starters_hash.size());
        for (auto const &imap: starters_hash) starter_nodes.push_back(imap.first);
    }

    std::vector<std::thread> threads;
    dbg_boss::size_type seqs_per_thread = starter_nodes.size() / n_threads;
    rem = starter_nodes.size() % n_threads;
    size_t amb_seqs=0, fm_seqs=0, nf_seqs=0, reb_seqs=0;

    std::cout<<"inferring the original sequences"<<std::endl;
    for (size_t i = 0; i < n_threads; i++) {

        start = i * seqs_per_thread;
        end = ((i + 1) * seqs_per_thread)-1;
        if (i == n_threads - 1) end += rem;
        threads.emplace_back(std::thread(build_sequences,
                                         start,
                                         end,
                                         std::ref(reb_seqs),
                                         std::ref(amb_seqs),
                                         std::ref(fm_seqs),
                                         std::ref(nf_seqs),
                                         std::ref(starter_nodes),
                                         std::ref(dbg_index),
                                         std::ref(fm_index),
                                         std::ref(output_file)));
    }
    for (auto &it : threads) it.join();
    output_file.close();
    std::cout<<reb_seqs<<" were rebuilt and "<<amb_seqs<<" were ambiguous, " << fm_seqs << " were fm_count 0, " << nf_seqs << " were not found"<<std::endl;
}

