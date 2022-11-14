#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/csa_wt.hpp>

const char termCharacter = '$';

void ParsePrint(std::string s, std::stringstream &ss) {
  std::string tmp_str;
  int count = 1;
  for (auto const &c : s) {
    if (c == '$') {
      if (tmp_str.size() == 0) continue;
      ss << ">id:" << count++ << '\n';
      ss << tmp_str << '\n';
      tmp_str = "";
    } else {
      tmp_str += c;
    }
  }
}

int main(int argc, char **argv) {
  // Usage: ./rebuild_st.out <index.cst> <output>"

  std::string index_path = std::string(argv[1]);
  std::string output_path = std::string(argv[2]) + ".fasta";

  sdsl::cst_sct3<> cst;
  sdsl::load_from_file(cst, index_path);

  std::stringstream ss;
  ParsePrint(sdsl::extract(cst.csa, 0, cst.csa.size() - 1), std::ref(ss));

  std::ofstream ofs(output_path);
  ofs << ss.str();
  ofs.close();

  return 0;
}