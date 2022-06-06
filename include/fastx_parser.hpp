#ifndef COL_BOSS_FASTX_PARSER_HPP
#define COL_BOSS_FASTX_PARSER_HPP

extern "C"{
#include "kseq.h"
#include <zlib.h>
}

#include <iostream>
#include <sdsl/sd_vector.hpp>

class fastx_parser {
    KSEQ_INIT(gzFile, gzread);
public:
    static int preproc_reads(std::string &input_file, sdsl::cache_config& config, size_t kmer_size);
};


#endif //COL_BOSS_FASTX_PARSER_HPP
