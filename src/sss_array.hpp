/*
  Part of SubseqSketch.
  File IO for binary sketching files and distance computation.
  By Ke @ Penn State
*/

#ifndef __SSS_ARRAY_H__
#define __SSS_ARRAY_H__

#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class sss_array
{
public:
    // Write a sketching matrix to file in binary format, dimension is
    // num_sketches x sketch_len, each row is the sketching of one sequence
    static void write_all(const Eigen::MatrixXi& sketches,
			  size_t num_sketches,
			  int sketch_len,
			  int max_val,
			  const std::string& sketch_file);

    // Load a binary file with an sketching matrix. The first three
    // values are assumbed to be num_sketches(size_t), sketch_len(int),
    // and max_val(int)
    static Eigen::MatrixXi load_all(size_t& num_sketches,
				    int& sketch_len,
				    int& max_val,
				    const std::string& sketch_file);
    
    // Write a single sketching array to file in binary format
    static void write(const int* sketch, int size, int max_val, std::ofstream& fout);
    
    // Read a binary file with an unknown number of sketchings, each
    // sketching array is preceded by its dimension, require all sketchings
    // in the file have the same dimension.
    static void load(std::vector<int*>& sketches, int& sketch_dim, int& num_tokens,
		     const std::string& sketch_file);

    // Compute row_normalized(sketch1) * column_normalized(transpose(sketch2))
    // and write the resulting matrix to dist_file in binary format.
    static void pairwise_cos_dist(const std::vector<int*>& sketch1,
				  const std::vector<int*>& sketch2,
				  int sketch_dim,
				  const std::string& dist_file);

    static void pairwise_cos_dist(const Eigen::MatrixXi& sketch1,
				  const Eigen::MatrixXi& sketch2_tran,
				  const std::string& dist_file);

    // Free each int array in sketches.
    static void free(std::vector<int*>& sketches);

    // Load a distance matrix saved by save_dist_matrix and output
    // to std::cout or to npy format
    static void load_dist_matrix(const std::string& dist_file, bool to_stdout);
    
private:
    static void save_dist_matrix(const Eigen::MatrixXd& dist,
				 const std::string& dist_file);
    static void show_dist_matrix(const Eigen::MatrixXd& dist);
    static void save_dist_matrix_to_npy(const Eigen::MatrixXd& dist,
					const std::string& dist_file);
};

#endif
