/*
  Part of SubseqEmbed.
  File IO for binary embedding files.
  By Ke @ Penn State
*/

#include "rssebd_array.hpp"
#include <string>
#include <Eigen/Dense>

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
	// throw std::runtime_error("Could not open the file: " + embed_file);
	std::cerr << "Error: could not open the file: "
		  << embed_file << std::endl;
	std::exit(1);
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
	    /*
	    throw std::runtime_error("Inconsistent embedding dimension found, #1: " +
				     std::to_string(embed_dim) + " #" +
				     std::to_string(embeds.size()+1) + ": " +
				     std::to_string(cur_dim));
	    */
	    std::cerr << "Error: inconsistent embedding dimension found, #1: "
		      << embed_dim << " #" << embeds.size() + 1
		      << ": " << cur_dim << std::endl;
	    std::exit(1);
	    
	}

	int cur_num_tokens;
	fin.read(reinterpret_cast<char*>(&cur_num_tokens), sizeof(cur_num_tokens));
	if(num_tokens < 0) num_tokens = cur_num_tokens;
	else if(num_tokens != cur_num_tokens)
	{
	    /*
	    throw std::runtime_error("Inconsistent max value found, #1: " +
				     std::to_string(num_tokens) + " #" +
				     std::to_string(embeds.size()+1) + ": " +
				     std::to_string(cur_num_tokens));
	    */
	    std::cerr << "Warning: inconsistent max value found, #1: "
		      << num_tokens << " #" << embeds.size() + 1
		      << ": " << cur_num_tokens << std::endl;
	}

	int* cur = new int[cur_dim];
	fin.read(reinterpret_cast<char*>(cur), sizeof *cur * embed_dim);
	embeds.push_back(cur);
    }

    fin.close();
}


void rssebd_array::pairwise_cos_dist(const std::vector<int*>& embed1,
				     const std::vector<int*>& embed2,
				     int embed_dim,
				     const std::string& dist_file)
{
    std::ofstream fout(dist_file);
    if(!fout)
    {
	// throw std::runtime_error("Could not write to the file: " + dist_file);
	std::cerr << "Error: could not write to the file: "
		  << dist_file << std::endl;
	std::exit(1);
    }

    std::size_t s1 = embed1.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m1(s1, embed_dim);
    for(int i = 0; i < s1; ++i)
    {
	for(int j = 0; j < embed_dim; ++j)
	{
	    m1(i, j) = embed1[i][j];
	}
    }
    m1.rowwise().normalize();

    std::size_t s2 = embed2.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m2(s2, embed_dim);
    for(int i = 0; i < s2; ++i)
    {
	for(int j = 0; j < embed_dim; ++j)
	{
	    m2(i, j) = embed2[i][j];
	}
    }
    m2.rowwise().normalize();
    m2.transposeInPlace();

    Eigen::MatrixXd dist(s1, s2);
    dist.noalias() = m1 * m2;
    
    Eigen::IOFormat format(Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n");

    fout << dist << std::endl;
    fout.close();
}


void rssebd_array::free(std::vector<int*>& embeds)
{
    for(int* x : embeds)
    {
	delete[] x;
    }
    embeds.clear();
}
