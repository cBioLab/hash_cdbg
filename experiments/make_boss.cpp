#include <iostream>
#include <build_index.hpp>

int main(int argc, char* argv[]) {
  std::string input_file = argv[1];
  size_t kmer_size = 30;

  std::cout << "k-mer: " << kmer_size << std::endl;

  cache_config config;
  build_index::build_boss(input_file, config, kmer_size);

  return 0;
}
