#include <chrono>
#include <fstream>
#include <iostream>
#include <mutex>
#include <regex>
#include <thread>
#include <vector>

struct fastq {
  std::string name;
  std::string read;
  std::string id;
  std::string score;
};

std::vector<std::vector<fastq>::iterator> div_itr(int n, std::vector<fastq>::iterator begin, std::vector<fastq>::iterator end) {
  std::vector<std::vector<fastq>::iterator> result;
  const int length = std::distance(begin, end);
  const int partial_length = (length + n - 1) / n;

  for (std::vector<fastq>::iterator i = begin; i < end; std::advance(i, partial_length)) {
    result.push_back(i);
  }
  result.push_back(end);

  return result;
}

void worker_remove(std::vector<fastq>& reads, std::vector<fastq>::iterator begin, std::vector<fastq>::iterator end, std::vector<fastq>& i_formated_reads) {
  static const std::regex regex(".*(n|N).*");

  for (auto itr = begin; itr != end; ++itr) {
    if (std::regex_search(itr->read, regex)) continue;
    i_formated_reads.push_back({itr->name, itr->read, itr->id, itr->score});
  }
}

int main(int argc, char const *argv[]) {
  /*
  usage: ./remove_n_read input.fastq output.fastq [thread_num]
  */

  std::chrono::system_clock::time_point  start, middle, end;
  double elapsed = 0;
  start = std::chrono::system_clock::now();

  const std::string infile = std::string(argv[1]);
  std::ifstream ifs(infile);
  if (!ifs) {
    std::cerr << "Cannot open " << infile << " !" << std::endl;
    return 0;
  }

  std::vector<fastq> reads;
  std::string input_name;
  std::string input_read;
  std::string input_id;
  std::string input_score;

  while (std::getline(ifs, input_name)) {
    std::getline(ifs, input_read);
    std::getline(ifs, input_id);
    std::getline(ifs, input_score);
    reads.push_back({input_name, input_read, input_id, input_score});
  }
  ifs.close();

  middle = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast < std::chrono::milliseconds > (middle - start).count();
  std::cout << "Middle. " << elapsed << "ms (" << elapsed / 60000 << "min)" << '\n';

  int thread_num = 8;
  if (argc > 3) thread_num = stoi(std::string(argv[3]));
  std::cout << "Using " << thread_num << " threads." << std::endl;
  std::vector<std::thread> ths(thread_num);
  std::vector<std::vector<fastq>::iterator> div_reads_itr = div_itr(thread_num, reads.begin(), reads.end());

  std::vector<std::vector<fastq>> formated_reads(thread_num);
  int i_itr = 0;
  for (auto& th : ths) {
    th = std::thread(worker_remove, std::ref(reads), div_reads_itr[i_itr], div_reads_itr[i_itr + 1], std::ref(formated_reads[i_itr]));
    i_itr++;
  }

  for (auto& th : ths) {
    th.join();
  }

  std::ofstream ofs((std::string(argv[2])));
  for (int i = 0; i < thread_num; ++i) {
    std::vector<fastq>::iterator reads_end = formated_reads[i].end();
    for (auto itr = formated_reads[i].begin(); itr != reads_end; ++itr) {
      ofs << itr->name << '\n';
      ofs << itr->read << '\n';
      ofs << itr->id << '\n';
      ofs << itr->score << '\n';
    }
  }
  ofs.close();

  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();
  std::cout << "Finish. " << elapsed << "ms (" << elapsed / 60000 << "min)" << '\n';

  return 0;
}
