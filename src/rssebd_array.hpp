/*
  Part of SubseqEmbed.
  File IO for binary embedding files.
  By Ke @ Penn State
*/

#ifndef __RSSEBD_ARRAY_H__
#define __RSSEBD_ARRAY_H__

#include <fstream>
#include <vector>

class rssebd_array
{
public:
    // Write a single embedding array to file in binary format
    static void write(const int* embed, int size, int max_val, std::ofstream& fout);
    
    // Read a binary file with an unknown number of embeddings, each
    // embedding array is preceded by its dimension, require all embeddings
    // in the file have the same dimension.
    static void load(std::vector<int*>& embeds, int& embed_dim, int& num_tokens,
		     const std::string& embed_file);
};

#endif
