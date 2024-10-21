# Overview

SubseqEmbed creates embeddings of sequences using random subsequences.

# Installation
Clone SubseqEmbed from [github](https://github.com/Shao-Group/SubseqEmbed.git)

## Compile SubseqEmbed
Use the following to compile SubseqEmbed:
```
cd SubseqEmbed
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```
This should create an executable `SubseqEmbed` within the `build/` directory.

# Usage
Check the usage information of SubseqEmbed and its subcommands:
```
build/SubseqEmbed [SUBCOMMAND] --help
```
