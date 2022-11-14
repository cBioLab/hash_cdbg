#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <mutex>
#include <map>

int main(int argc, char const *argv[]) {
  if (argc <= 2) {
    std::cerr << "usage: ./check_rebuild_read original.fa rebuild.fa" << std::endl;
    return 0;
  }

  std::ifstream ifs_original((std::string(argv[1])));
  if (!ifs_original) {
    std::cerr << "Cannot open " << std::string(argv[1]) << " !" << std::endl;
    return 0;
  }
  std::ifstream ifs_rebuild((std::string(argv[2])));
  if (!ifs_rebuild) {
    std::cerr << "Cannot open " << std::string(argv[2]) << " !" << std::endl;
    return 0;
  }

  std::map<std::string, long long> original_read;
  std::string name, read;
  long long count = 0;
  while(std::getline(ifs_original, name)){
    if ((name.at(0) == '@' || name.at(0) == '>') && std::getline(ifs_original, read)) {
      std::string rc_read;
      for (int i = read.size() - 1; i >= 0; --i) {
        switch (read[i]) {
          case 'A':
          case 'a':
            rc_read += 'T';
            break;
          case 'T':
          case 't':
            rc_read += 'A';
            break;
          case 'G':
          case 'g':
            rc_read += 'C';
            break;
          case 'C':
          case 'c':
            rc_read += 'G';
            break;
          default:
            rc_read += 'N';
            break;
        }
      }
      original_read[read]++;
      original_read[rc_read]++;
      count += 2;
    }
  }
  ifs_original.close();
  std::cout << "map size: " << original_read.size() << ", original read size: " << count << std::endl;

  std::map<std::string, long long> rebuild_read;
  count = 0;
  while(std::getline(ifs_rebuild, name)){
    if ((name.at(0) == '@' || name.at(0) == '>') && std::getline(ifs_rebuild, read)) {
      if (original_read.find(read) == original_read.end()) {
        rebuild_read[read]++;
      } else {
        original_read[read]--;
      }
      count++;
    }
  }
  ifs_rebuild.close();
  std::cout << "rebuild read size: " << count << std::endl;

  long long count_nr = 0, count_no = 0;
  std::cout << "### Not rebuild reads ###" << std::endl;
  for (auto itr = original_read.begin(); itr != original_read.end(); ++itr) {
    if (itr->second > 0) {
      std::cout << itr->first << "\t" << itr->second << std::endl;
      count_nr += itr->second;
    } else if (itr->second < 0) {
      rebuild_read[itr->first] -= itr->second;
    }
  }

  std::cout << "### Not in original read ###" << std::endl;
  for (auto itr = rebuild_read.begin(); itr != rebuild_read.end(); ++itr) {
    if (itr->second > 0) {
      std::cout << itr->first << "\t" << itr->second << std::endl;
      count_no += itr->second;
    }
  }

  std::cout << "Not rebuild: " << count_nr << " reads, Not in original: " << count_no << " reads" << std::endl;

  return 0;
}