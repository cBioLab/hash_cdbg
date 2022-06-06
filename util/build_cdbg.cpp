#include <iostream>
#include <hash_cdbg/boss.hpp>

int main(int argc, char* argv[]) {
  std::string input_file = "data/example.fastq";
  size_t kmer_size = 30;
  size_t n_threads = 1;

  dbg_boss dbg_index(input_file, kmer_size, n_threads);
  store_to_file(dbg_index, "example.cdbg");

  return 0;
}
