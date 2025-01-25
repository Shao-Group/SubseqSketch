# Overview

SubseqSketch creates sketchings of sequences using random subsequences.
It exhibits strong correlation with the edit distance and outperforms other sketching methods in tasks such as nearest neighbor search and phylogeny reconstruction.
See [SubseqSketch-test](https://github.com/Shao-Group/SubseqSketch-test.git) for scripts that can be used to reproduce the evaluation results in the SubseqSketch paper.

SubseqSketch uses [Eigen](https://gitlab.com/libeigen/eigen) for matrix multiplication, the needed header files are included in this repository.
Optionally, SubseqSketch uses [OpenMP](https://www.openmp.org/resources/openmp-compilers-tools/) for parallel computing.

# Installation
Clone SubseqSketch from [github](https://github.com/Shao-Group/SubseqSketch.git)

## Compile SubseqSketch
SubseqSketch requires [CMake](https://cmake.org/download/) to build.
Use the following to compile SubseqSketch:
```
cd SubseqSketch
mkdir build && cd build
cmake ..
make -j
```
This should create an executable `SubseqSketch` within the `build/` directory.

# Usage
A general pipeline for sketching with SubseqSketch is provided below.
The detailed parameter information of SubseqSketch and its subcommands can be checked by
```
build/SubseqSketch [SUBCOMMAND] --help
```
1. Generate a list of $n$ random testing subsequences with token size $t$, each subsequences contain $l$ tokens. They will be stored in a plain text file `subsequences.txt` by default.
   ```
   build/SubseqSketch init -a alphabets/DNA -n 128 -t 3 -l 15
   ```
2. Use the random testing subsequences to generate SubseqSketches for sequences in fasta format:
   ```
   build/SubseqSketch sketch -s subsequences.txt input1.fa input2.fa ...
   ```
   For each input file the sketches will be stored in a binary file such as `input1.n128.l15.t3.sss` which can be viewed (or redirect to a plain text file) by the `info` subcommand:
   ```
   build/SubseqSketch info input1.n128.l15.t3.sss | less
   ```
3. Compute the all-vs-all sketch distances between two sketches:
   ```
   build/SubseqSketch dist -o input1-vs-input2.sss-dist input1.n128.l15.t3.sss input2.n128.l15.t3.sss
   ```
   If `input1.fa` has $s_1$ sequences and `input2.fa` has $s_2$ sequences, then the result is a $s_1\times s_2$ matrix $M$ where $M_{i,j}$ is the cosine distance between (the sketches of) the $i$-th sequence in `input1.fa` and the $j$-th sequence in `input2.fa`.
   The matrix is stored in binary format, which can be viewed by:
   ```
   build/SubseqSketch show input1-vs-input2.sss-dist | less
   ```
   or transformed into a npy format for downstream tasks using `numpy`:
   ```
   build/SubseqSketch show -p input1-vs-input2.sss-dist
   ```
