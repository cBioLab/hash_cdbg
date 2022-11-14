#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/csa_wt.hpp>

const char termCharacter = '$';

char getRCBase(char base) {
  switch(base) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    default : return 'N';
  }
}

std::string getRCRead(const std::string& read) {
  int readlen = read.size();
  std::string result;
  for(int i = 0; i < readlen; i++) {
    result += getRCBase(read[i]);
  }
  reverse(result.begin(), result.end());
  return result;
}

bool isValidRead(const std::string &read) {
  const std::string alphabet = "ACGT";
  int readlen = read.size();
  for(int i = 0; i < readlen; i++){
    if(std::find(alphabet.begin(), alphabet.end(), read[i]) == alphabet.end()) return false;
  }
  return true;
}

int main(int argc, char **argv) {
  // Usage: ./make_st.out <read.fastq> <index-file-name> <thread-num>"

  std::string readFilePath = std::string(argv[1]);
  std::string indexFilePath = std::string(argv[2]);
  std::string tmpFilePath = indexFilePath + "_read_text.tmp";
  std::ifstream ifsRead(readFilePath);
  std::ofstream ofsTmp(tmpFilePath);

  ofsTmp << termCharacter;
  std::string name, seq, blankline, qual;
  while(getline(ifsRead, name)) {
    getline(ifsRead, seq);
    getline(ifsRead, blankline);
    getline(ifsRead, qual);

    if(!isValidRead(seq)) continue;

    ofsTmp << seq << termCharacter;
    std::string rcseq = getRCRead(seq);
    ofsTmp << rcseq << termCharacter;
  }
  ifsRead.close();
  ofsTmp.close();

  std::string CSTFilePath = indexFilePath + ".cst";
  // std::string CSAFilePath = indexFilePath + ".csa";

  sdsl::cst_sct3<> cst;
  // sdsl::csa_wt<> csa;

  sdsl::memory_monitor::start();
  auto start = std::chrono::high_resolution_clock::now();

  sdsl::cache_config config;
  construct(cst, tmpFilePath, config, 1);
  store_to_file(cst, CSTFilePath);

  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Finish CST. " << std::chrono::duration_cast<std::chrono::seconds>(stop-start).count() << " seconds." << std::endl;
  std::cout << "Peak memory by SDSL: " << sdsl::memory_monitor::peak() / 1024 << " KB" << std::endl;

  stop = std::chrono::high_resolution_clock::now();
  std::cout << "Finish. " << std::chrono::duration_cast<std::chrono::seconds>(stop-start).count() << " seconds." << std::endl;

  sdsl::memory_monitor::stop();
  std::cout << "Peak memory by SDSL: " << sdsl::memory_monitor::peak() / 1024 << " KB" << std::endl;

  return 0;
}