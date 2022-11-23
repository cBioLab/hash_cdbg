extern "C"{
#include "boss/kseq.h"
#include <zlib.h>
};

#include <sdsl/sd_vector.hpp>
#include <sdsl/suffix_arrays.hpp>
#include "dna_alphabet.hpp"


KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " INPUT_FASTQ OUTPUT_PREFIX" << std::endl;
        return 1;
    }

    cache_config config(false, ".", "tmp");
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;

    uint8_t int_symbol;
    string output_file = cache_file_name(conf::KEY_TEXT, config);

    int_vector_buffer<8> int_text(output_file, std::ios::out, 1000000);

    gzFile fp;
    kseq_t *seq;
    int l;
    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);
    size_t i = 0;

    int_text[i] = dna_alphabet::char2comp['$'];
    i++;
    size_t cont=0;
    while ((l = kseq_read(seq)) >= 0) {
        for (size_t j = 0; j < l; j++) {
            int_symbol = dna_alphabet::char2comp[seq->seq.s[j]];
            int_text[i] = int_symbol;
            i++;
        }
        int_text[i] = dna_alphabet::char2comp['$'];
        i++;
        cont++;
    }
    kseq_destroy(seq);
    gzclose(fp);
    int_text.close();

    construct(fm_index, cache_file_name(conf::KEY_TEXT, config));
    store_to_file(fm_index, string(argv[2]) + ".fm_index");
}
