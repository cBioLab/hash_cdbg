# hash_cdbg
A C++ library for indexing genome sequencing datasets by using colored de Bruijn Graph, hash functions and Bloom Filter.
The implementation is based on [this](https://bitbucket.org/DiegoDiazDominguez/colored_bos/src/master/) library by Diego Diaz Dominguez et al.

## Requirements
This tool requires:
- [CMake](https://cmake.org/) 3.15 or higher
- [SDSL-lite](https://github.com/simongog/sdsl-lite) library

## Support
- Ubuntu 18.04

## Installation

First, download the library and move to library's root directory.

```bash
git clone git@github.com:cBioLab/hash_cdbg.git
cd hash_cdbg
```

Then, prepare for compilation.

```bash
mkdir build && cd build
cmake ..
```

If you want to specify the directory in which to install this library, you can use:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX={your_install_path}/hash_cdbg
```

Finally, compile and install the library.

```bash
make & make install
```

## Getting Started
To use this library quickly, look in the util directory.
```build_cdbg.cpp``` is a code that builds an index, the detail of which is as follow:

```cpp
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
```

To compile and execute this code, do the following:

```bash
g++ -o build_cdbg.out ./util/build_cdbg.cpp -I {your_install_path}/include -L {your_install_path}/lib -lhash_cdbg -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -lz -std=c++17 -O3
./build_cdbg.out
```

The resulting ```example.cdbg``` is the index file.
To rebuild the original sequences from this index, do the following using ```build_fm_index.cpp``` and ```rebuild_seqs.cpp```:

```bash
g++ -o build_fm_index.out ./util/build_fm_index.cpp -I {your_install_path}/include -L {your_install_path}/lib -lhash_cdbg -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -lz
./build_fm_index.out data/example.fastq example
g++ -o rebuild_seqs.out ./util/rebuild_seqs.cpp -I {your_install_path}/include -L {your_install_path}/lib -lhash_cdbg -lsdsl -ldivsufsort -ldivsufsort64 -lpthread -lz -std=c++17 -O3
./rebuild_seqs.out example.cdbg example.fm_index 1 example.re
```

The resulting ```example.re.fasta``` is a fasta file that contains the ```example.fastq``` sequences and it's reverse complements rebuilt.

## Reproduction of Our Experiments
If you want to reproduce our experiments, see [experiments README](experiments/README.md).
