#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char const *argv[]) {
  if (argc <= 2) {
    cerr << "Using: ./Fastq2Fasta <read.fastq> <output.fasta>" << endl;
    return 0;
  }

  ifstream ifs((string(argv[1])));
  if (!ifs) {
    cerr << "Cannot open " << string(argv[1]) << " !" << endl;
    return 0;
  }
  ofstream ofs((string(argv[2])));
  string input;

  while(getline(ifs, input)) {
    if (input[0] == '>' || input[0] == '@') {
      ofs << '>' << input.erase(0, 1) << '\n';
      if (getline(ifs, input)) {
        ofs << input << '\n';
      } else {
        ofs <<'\n';
      }
      getline(ifs, input);
      getline(ifs, input);
    }
  }

  ifs.close();
  ofs.close();

  return 0;
}