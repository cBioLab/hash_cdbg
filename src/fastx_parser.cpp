#include "dna_alphabet.hpp"
#include "fastx_parser.hpp"

using namespace std;
using namespace sdsl;

int fastx_parser::preproc_reads(string &input_file, cache_config& config, size_t kmer_size) {

    //TODO consider complementary reverse and handle N values with $ symbols
    gzFile fp;
    kseq_t *seq;
    size_t len=0, n_seq=0;
    string output_file = cache_file_name(conf::KEY_TEXT, config);
    string dollar_bv_file = cache_file_name("dollar_bv", config);

    int_vector_buffer<8> int_text(output_file, std::ios::out, 1000000, 3);
    int_vector_buffer<1> dollar_bv(dollar_bv_file, std::ios::out, 100000, 1);

    fp = gzopen(input_file.c_str(), "r");
    seq = kseq_init(fp);
    size_t i=0;

    ofstream out(output_file);
    int_text[i] = dna_alphabet::char2comp['$'];
    dollar_bv[i] = true;
    i++;

    int_vector_buffer<8>::value_type int_symbol;

    while(kseq_read(seq)>=0){
        if(n_seq==0){
            //TODO check if all sequences have the same length (?)
            len = strlen(seq->seq.s);
            //TODO report number of sequences read (?)
            n_seq++;
        }

        for(size_t j=len;j --> 0;){
            // I have to handle N symbols
            if(seq->seq.s[j] == 'N') {
                int_symbol = dna_alphabet::char2comp['$'];
                dollar_bv[i] = true;
            }else{
                int_symbol = dna_alphabet::char2comp[seq->seq.s[j]];
                dollar_bv[i] = false;
            }
            int_text[i] = int_symbol;
            i++;
        }
        //add dollar symbols
        for(size_t j=0;j<kmer_size;j++){
            int_text[i] = dna_alphabet::char2comp['$'];
            dollar_bv[i] = true;
            i++;
        }
        //add complementary reverse sequences
        for(size_t j=0;j<len; j++){
            // I have to handle N symbols
            if(seq->seq.s[j] == 'N') {
                int_symbol = dna_alphabet::char2comp['$'];
                dollar_bv[i] = true;
            }else{
                int_symbol = dna_alphabet::comp2rev[dna_alphabet::char2comp[seq->seq.s[j]]];
                dollar_bv[i] = false;
            }
            int_text[i] = int_symbol;
            i++;
        }
        //add dollar symbols again
        for(size_t j=0;j<kmer_size;j++){
            int_text[i] = dna_alphabet::char2comp['$'];
            dollar_bv[i] = true;
            i++;
        }
    }

    //least symbol
    int_text[i]=0;
    dollar_bv[i] = false;

    int_text.close();
    dollar_bv.close();

    bit_vector marked_dollars;
    load_from_file(marked_dollars, dollar_bv_file);
    sd_vector<> sp_dollar_bv(marked_dollars);
    store_to_file(sp_dollar_bv, cache_file_name("sp_dollar_bv", config));

    if(remove(dollar_bv_file)!=0){
       perror("Error removing bit vector file");
    }

    register_cache_file(conf::KEY_TEXT, config);
    register_cache_file("sp_dollar_bv", config);

    kseq_destroy(seq);
    gzclose(fp);
    return EXIT_SUCCESS;
}

/*int fastx_parser::init_map_file(char **map_seq, off_t file_size, int *fd, string output_file) {

    off_t result;
    *fd = open(output_file.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);

    if (*fd == -1) {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    result = lseek(*fd, file_size - 1, SEEK_SET);
    if (result == -1) {
        close(*fd);
        perror("Error calling lseek() to 'stretch' the file");
        exit(EXIT_FAILURE);
    }

    result = write(*fd, "", 1);
    if (result != 1) {
        close(*fd);
        perror("Error writing last byte of the file");
        exit(EXIT_FAILURE);
    }

    *map_seq = static_cast<char *>(mmap(nullptr, static_cast<size_t>(file_size), PROT_READ | PROT_WRITE,
                                       MAP_SHARED, *fd, 0));
    if (*map_seq == MAP_FAILED) {
        close(*fd);
        perror("Error mapping the file");
        exit(EXIT_FAILURE);
    }

    return 0;
}*/

