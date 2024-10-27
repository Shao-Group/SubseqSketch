/*
  Part of SubseqEmbed.
  File IO for binary embedding files and distance computation.
  By Ke @ Penn State
*/

#ifndef __RSSEBD_ARRAY_H__
#define __RSSEBD_ARRAY_H__

#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

class rssebd_array
{
public:
    // Write a embedding matrix to file in binary format, dimension is
    // num_embeds x embed_len, each row is the embedding of one sequence
    static void write_all(const Eigen::MatrixXi& embeds, size_t num_embeds,
			  int embed_len, int max_val, const std::string& embed_file);

    // Load a binary file with an embedding matrix. The first three
    // values are assumbed to be num_embeds(size_t), embed_len(int),
    // and max_val(int)
    static Eigen::MatrixXd load_all(size_t& num_embeds,
				    int& embed_len,
				    int& max_val,
				    bool normalize,
				    bool transpose,
				    const std::string& embed_file);
    
    // Write a single embedding array to file in binary format
    static void write(const int* embed, int size, int max_val, std::ofstream& fout);
    
    // Read a binary file with an unknown number of embeddings, each
    // embedding array is preceded by its dimension, require all embeddings
    // in the file have the same dimension.
    static void load(std::vector<int*>& embeds, int& embed_dim, int& num_tokens,
		     const std::string& embed_file);

    // Compute row_normalized(embed1) * column_normalized(transpose(embed2))
    // and write the resulting matrix to dist_file in binary format.
    static void pairwise_cos_dist(const std::vector<int*>& embed1,
				  const std::vector<int*>& embed2,
				  int embed_dim,
				  const std::string& dist_file);

    // Assume embed2 has been transposed, see rssebd_array::load_all.
    static void pairwise_cos_dist(const Eigen::MatrixXd& embed1,
				  const Eigen::MatrixXd& embed2_tran,
				  const std::string& dist_file);

    // Assume embed2 has NOT been transposed, see rssebd_array::load_all.
    static void pairwise_max_likelyhood_dist(const Eigen::MatrixXd& embed1,
					     const Eigen::MatrixXd& embed2,
					     const std::string& dist_file);

    // Free each int array in embeds.
    static void free(std::vector<int*>& embeds);

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
