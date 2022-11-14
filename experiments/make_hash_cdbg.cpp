#include <chrono>
#include <iostream>

#include <hash_cdbg/boss.hpp>

int main(int argc, char const *argv[]){
  std::cout << "command: $ ";
  for (int i = 0; i < argc; ++i) {
    std::cout << argv[i] << " ";
  }
  std::cout << std::endl;

  std::string file = argv[1];  // it can also be a fasta
  std::string out = argv[2];
  size_t kmer_size = 30;
  size_t num_thread = 16;

  std::chrono::system_clock::time_point  start, end;
  double elapsed = 0;
  start = std::chrono::system_clock::now();

  sdsl::memory_monitor::start();

  dbg_boss dbg_index(file, kmer_size, num_thread);
  store_to_file(dbg_index, out);

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();
  std::cout << "Finish. " << elapsed << "ms (" << elapsed / 60000 << "min)" << '\n';

  sdsl::memory_monitor::stop();
  std::cout << "Peak memory by SDSL: " << sdsl::memory_monitor::peak() / 1024 << " KB" << std::endl;

  return 0;
}
