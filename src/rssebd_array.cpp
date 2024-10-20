/*
  Part of SubseqEmbed.
  File IO for binary embedding files.
  By Ke @ Penn State
*/

#include "rssebd_array.hpp"
#include <string>

void rssebd_array::write(const int* embed, int size, int max_val,
			 std::ofstream& fout)
{
    fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
    fout.write(reinterpret_cast<const char*>(&max_val), sizeof(max_val));
    fout.write(reinterpret_cast<const char*>(embed), sizeof *embed * size);
}

void rssebd_array::load(std::vector<int*>& embeds, int& embed_dim,
			int& num_tokens, const std::string& embed_file)
{
    std::ifstream fin(embed_file, std::ios::binary);

    if(!fin)
    {
	throw std::runtime_error("Could not open the file: " + embed_file);
    }

    embed_dim = -1;
    num_tokens = -1;
    
    while(true)
    {
	int cur_dim;
	fin.read(reinterpret_cast<char*>(&cur_dim), sizeof(cur_dim));
	if(fin.eof()) break;
	
	if(embed_dim < 0) embed_dim = cur_dim;
	else if(embed_dim != cur_dim)
	{
	    throw std::runtime_error("Inconsistent embedding dimension found, #1: " +
				     std::to_string(embed_dim) + " #" +
				     std::to_string(embeds.size()+1) + ": " +
				     std::to_string(cur_dim));
	}

	int cur_num_tokens;
	fin.read(reinterpret_cast<char*>(&cur_num_tokens), sizeof(cur_num_tokens));
	if(num_tokens < 0) num_tokens = cur_num_tokens;
	else if(num_tokens != cur_num_tokens)
	{
	    throw std::runtime_error("Inconsistent max value found, #1: " +
				     std::to_string(num_tokens) + " #" +
				     std::to_string(embeds.size()+1) + ": " +
				     std::to_string(cur_num_tokens));
	}

	int* cur = new int[cur_dim];
	fin.read(reinterpret_cast<char*>(cur), sizeof *cur * embed_dim);
	embeds.push_back(cur);
    }

    fin.close();
}
