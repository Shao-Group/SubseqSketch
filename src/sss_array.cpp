/*
  Part of SubseqSketch.
  File IO for binary sketching files.
  By Ke @ Penn State
*/

#include "sss_array.hpp"
#include <string>

void sss_array::write_all(const Eigen::MatrixXi& sketches, size_t num_sketches,
			  int sketch_len, int max_val,
			  const std::string& sketch_file)
{
    std::ofstream fout(sketch_file, std::ios::binary);
    if(!fout)
    {
	std::cerr << "Error: could not open the file: "
		  << sketch_file << std::endl;
	std::exit(1);
    }
    
    fout.write(reinterpret_cast<const char*>(&num_sketches), sizeof(num_sketches));
    fout.write(reinterpret_cast<const char*>(&sketch_len), sizeof(sketch_len));
    fout.write(reinterpret_cast<const char*>(&max_val), sizeof(max_val));

    fout.write(reinterpret_cast<const char*>(sketches.data()), sizeof(int) * num_sketches * sketch_len);
    fout.close();
}

Eigen::MatrixXi
sss_array::load_all(size_t& num_sketches,
		    int& sketch_len,
		    int& max_val,
		    const std::string& sketch_file)
{
    std::ifstream fin(sketch_file, std::ios::binary);

    if(!fin)
    {
	// throw std::runtime_error("Could not open the file: " + sketch_file);
	std::cerr << "Error: could not open the file: "
		  << sketch_file << std::endl;
	std::exit(1);
    }
    
    fin.read(reinterpret_cast<char*>(&num_sketches), sizeof(num_sketches));
    fin.read(reinterpret_cast<char*>(&sketch_len), sizeof(sketch_len));
    fin.read(reinterpret_cast<char*>(&max_val), sizeof(max_val));

    Eigen::MatrixXi sketches(num_sketches, sketch_len);

    fin.read(reinterpret_cast<char*>(sketches.data()), sizeof(int) * num_sketches * sketch_len);
    fin.close();

    return sketches;
}


void sss_array::write(const int* sketch, int size, int max_val,
		      std::ofstream& fout)
{
    fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
    fout.write(reinterpret_cast<const char*>(&max_val), sizeof(max_val));
    fout.write(reinterpret_cast<const char*>(sketch), sizeof *sketch * size);
}

void sss_array::load(std::vector<int*>& sketches, int& sketch_dim,
		     int& num_tokens, const std::string& sketch_file)
{
    std::ifstream fin(sketch_file, std::ios::binary);

    if(!fin)
    {
	// throw std::runtime_error("Could not open the file: " + sketch_file);
	std::cerr << "Error: could not open the file: "
		  << sketch_file << std::endl;
	std::exit(1);
    }

    sketch_dim = -1;
    num_tokens = -1;
    
    while(true)
    {
	int cur_dim;
	fin.read(reinterpret_cast<char*>(&cur_dim), sizeof(cur_dim));
	if(fin.eof()) break;
	
	if(sketch_dim < 0) sketch_dim = cur_dim;
	else if(sketch_dim != cur_dim)
	{
	    /*
	      throw std::runtime_error("Inconsistent sketching dimension found, #1: " +
	      std::to_string(sketch_dim) + " #" +
	      std::to_string(sketches.size()+1) + ": " +
	      std::to_string(cur_dim));
	    */
	    std::cerr << "Error: inconsistent sketching dimension found, #1: "
		      << sketch_dim << " #" << sketches.size() + 1
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
	      std::to_string(sketches.size()+1) + ": " +
	      std::to_string(cur_num_tokens));
	    */
	    std::cerr << "Warning: inconsistent max value found, #1: "
		      << num_tokens << " #" << sketches.size() + 1
		      << ": " << cur_num_tokens << std::endl;
	}

	int* cur = new int[cur_dim];
	fin.read(reinterpret_cast<char*>(cur), sizeof *cur * sketch_dim);
	sketches.push_back(cur);
    }

    fin.close();
}


void sss_array::pairwise_cos_dist(const std::vector<int*>& sketch1,
				  const std::vector<int*>& sketch2,
				  int sketch_dim,
				  const std::string& dist_file)
{
    std::size_t s1 = sketch1.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m1(s1, sketch_dim);
    for(int i = 0; i < s1; ++i)
    {
	for(int j = 0; j < sketch_dim; ++j)
	{
	    m1(i, j) = sketch1[i][j];
	}
    }
    m1.rowwise().normalize();

    std::size_t s2 = sketch2.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> m2(s2, sketch_dim);
    for(int i = 0; i < s2; ++i)
    {
	for(int j = 0; j < sketch_dim; ++j)
	{
	    m2(i, j) = sketch2[i][j];
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


void sss_array::pairwise_cos_dist(const Eigen::MatrixXi& sketch1,
				  const Eigen::MatrixXi& sketch2,
				  const std::string& dist_file)
{
    Eigen::MatrixXd normalized1 = sketch1.cast<double>();
    normalized1.rowwise().normalize();
    Eigen::MatrixXd normalized2 = sketch2.cast<double>();
    normalized2.rowwise().normalize();
    normalized2.transposeInPlace();
    
    Eigen::MatrixXd dist(sketch1.rows(), sketch2.rows());
    dist.noalias() = normalized1 * normalized2;
    dist.array() = 1 - dist.array();
    double zero_threshold = 1e-8;
    dist = (dist.array() < zero_threshold).select(0.0f, dist);

    save_dist_matrix(dist, dist_file);
}


void sss_array::save_dist_matrix(const Eigen::MatrixXd& dist,
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

void sss_array::load_dist_matrix(const std::string& dist_file, bool to_stdout)
{
    std::ifstream fin(dist_file, std::ios::binary);

    if(!fin)
    {
	// throw std::runtime_error("Could not open the file: " + sketch_file);
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


void sss_array::free(std::vector<int*>& sketches)
{
    for(int* x : sketches)
    {
	delete[] x;
    }
    sketches.clear();
}

void sss_array::show_dist_matrix(const Eigen::MatrixXd& dist)
{
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

    std::cout << dist.format(fmt) << std::endl;
}

void sss_array::save_dist_matrix_to_npy(const Eigen::MatrixXd& dist,
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
