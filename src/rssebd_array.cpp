/*
  Part of SubseqEmbed.
  File IO for binary embedding files.
  By Ke @ Penn State
*/

#include "rssebd_array.hpp"
#include <string>

void rssebd_array::write_all(const Eigen::MatrixXi& embeds, size_t num_embeds,
			     int embed_len, int max_val,
			     const std::string& embed_file)
{
    std::ofstream fout(embed_file, std::ios::binary);
    if(!fout)
    {
	std::cerr << "Error: could not open the file: "
		  << embed_file << std::endl;
	std::exit(1);
    }
    
    fout.write(reinterpret_cast<const char*>(&num_embeds), sizeof(num_embeds));
    fout.write(reinterpret_cast<const char*>(&embed_len), sizeof(embed_len));
    fout.write(reinterpret_cast<const char*>(&max_val), sizeof(max_val));

    fout.write(reinterpret_cast<const char*>(embeds.data()), sizeof(int) * num_embeds * embed_len);
    fout.close();
}


Eigen::MatrixXd rssebd_array::load_all(size_t& num_embeds, int& embed_len,
				       int& max_val, bool transpose,
				       const std::string& embed_file)
{
    std::ifstream fin(embed_file, std::ios::binary);

    if(!fin)
    {
	// throw std::runtime_error("Could not open the file: " + embed_file);
	std::cerr << "Error: could not open the file: "
		  << embed_file << std::endl;
	std::exit(1);
    }
    
    fin.read(reinterpret_cast<char*>(&num_embeds), sizeof(num_embeds));
    fin.read(reinterpret_cast<char*>(&embed_len), sizeof(embed_len));
    fin.read(reinterpret_cast<char*>(&max_val), sizeof(max_val));

    Eigen::MatrixXi embeds(num_embeds, embed_len);

    fin.read(reinterpret_cast<char*>(embeds.data()), sizeof(int) * num_embeds * embed_len);
    fin.close();

    Eigen::MatrixXd normalized = embeds.cast<double>();
    normalized.rowwise().normalize();

    if(transpose)
    {
	normalized.transposeInPlace();
    }

    return normalized;
}


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
    dist.array() = 1 - dist.array();
    double zero_threshold = 1e-8;
    dist = (dist.array() < zero_threshold).select(0.0f, dist);

    save_dist_matrix(dist, dist_file);
}


void rssebd_array::pairwise_cos_dist(const Eigen::MatrixXd& embed1,
				     const Eigen::MatrixXd& embed2_tran,
				     const std::string& dist_file)
{
    Eigen::MatrixXd dist(embed1.rows(), embed2_tran.cols());
    dist.noalias() = embed1 * embed2_tran;
    dist.array() = 1 - dist.array();
    double zero_threshold = 1e-8;
    dist = (dist.array() < zero_threshold).select(0.0f, dist);

    save_dist_matrix(dist, dist_file);
}


void rssebd_array::save_dist_matrix(const Eigen::MatrixXd& dist,
				    const std::string& dist_file)
{
    std::ofstream fout(dist_file, std::ios::binary);
    if(!fout)
    {
	// throw std::runtime_error("Could not write to the file: " + dist_file);
	std::cerr << "Error: could not write to the file: "
		  << dist_file << std::endl;
	std::exit(1);
    }

    int rows = dist.rows();
    int cols = dist.cols();
    fout.write(reinterpret_cast<char*>(&rows), sizeof(rows));
    fout.write(reinterpret_cast<char*>(&cols), sizeof(cols));

    fout.write(reinterpret_cast<const char*>(dist.data()), sizeof(double) * dist.size());
    fout.close();
}

void rssebd_array::load_dist_matrix(const std::string& dist_file, bool to_stdout)
{
    std::ifstream fin(dist_file, std::ios::binary);

    if(!fin)
    {
	// throw std::runtime_error("Could not open the file: " + embed_file);
	std::cerr << "Error: could not open the file: "
		  << dist_file << std::endl;
	std::exit(1);
    }

    int rows, cols;
    fin.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    fin.read(reinterpret_cast<char*>(&cols), sizeof(cols));

    Eigen::MatrixXd dist(rows, cols);

    fin.read(reinterpret_cast<char*>(dist.data()), sizeof(double) * dist.size());
    fin.close();

    std::cout << "Loaded " << rows << "x" << cols << " distance matrix" << std::endl;

    if(to_stdout)
    {
	show_dist_matrix(dist);
    }
    else
    {
	save_dist_matrix_to_npy(dist, dist_file + ".npy");
    }
}


void rssebd_array::free(std::vector<int*>& embeds)
{
    for(int* x : embeds)
    {
	delete[] x;
    }
    embeds.clear();
}

void rssebd_array::show_dist_matrix(const Eigen::MatrixXd& dist)
{
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

    std::cout << dist.format(fmt) << std::endl;
}

void rssebd_array::save_dist_matrix_to_npy(const Eigen::MatrixXd& dist,
					   const std::string& dist_file)
{
    std::ofstream fout(dist_file, std::ios::binary);
    if(!fout)
    {
	std::cerr << "Error: could not write to the file: "
		  << dist_file << std::endl;
	std::exit(1);
    }

    /*
      The first 6 bytes are a magic string: exactly \x93NUMPY.

      The next 1 byte is an unsigned byte: the major version number of the file format, e.g. \x01.

      The next 1 byte is an unsigned byte: the minor version number of the file format, e.g. \x00. Note: the version of the file format is not tied to the version of the numpy package.

      The next 2 bytes form a little-endian unsigned short int: the length of the header data HEADER_LEN.

      The next HEADER_LEN bytes form the header data describing the array’s format. It is an ASCII string which contains a Python literal expression of a dictionary. It is terminated by a newline (\n) and padded with spaces (\x20) to make the total of len(magic string) + 2 + len(length) + HEADER_LEN be evenly divisible by 64 for alignment purposes.
    */
    
    // NPY magic number and version
    const char npy_magic[] = { '\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0};
    fout.write(npy_magic, sizeof(npy_magic));

    // Create the header
    std::string header = "{'descr': '<f8', 'fortran_order': True, 'shape': (" 
	+ std::to_string(dist.rows()) + ", "
	+ std::to_string(dist.cols()) + "), }";

    int len = sizeof(npy_magic) + 2 + header.size();
    int padding_len =  (64 - (header.size() % 64)) % 64;
    std::string padding(padding_len, '\x20');

    // Write header length
    uint16_t header_len = static_cast<uint16_t>(header.size() + padding_len + 1);
    uint8_t header_len_le[2] = {static_cast<uint8_t>(header_len & 0xff),
				static_cast<uint8_t>((header_len >> 8) & 0xff)};
    fout.write(reinterpret_cast<char*>(header_len_le), 2);

    // Write header with padding
    fout << header << padding << '\n';

    // Write the data
    fout.write(reinterpret_cast<const char*>(dist.data()), dist.size() * sizeof(double));
    fout.close();

    std::cout << "Distance matrix wrote to the file: " << dist_file << std::endl;
}
