/*
  Part of SubseqEmbed.
  Edit distance embedding by random subsequences.
  By Ke @ Penn State
*/

#include <iostream>
#include <string>
#include <vector>

#include "fasta_reader.hpp"
#include "subsequences.hpp"
// #include "tokenized_sequence.hpp"
#include "rssebd_array.hpp"
#include "CLI11.hpp"

#include <omp.h>


void gen_random_subsequences(const int subseq_len,
			     const int token_len,
			     const int num_subseqs,
			     const std::string& alphabet_file,
			     const std::vector<std::string>& input_files,
			     const std::string& subseq_file);

void compute_embeddings(const std::string& subseq_file,
			const std::vector<std::string>& input_files);

void compute_distances(const std::string& embed_file1,
		       const std::string& embed_file2,
		       const std::string& dist_file);

void show_embeddings(const std::string& embed_file);

void show_distances(const std::string& dist_file, bool to_stdout);

int main(int argc, char** argv)
{
    CLI::App app("SubseqEmbed - Edit distance embedding by random subsequences");
    app.require_subcommand(1);
    app.get_formatter()->column_width(20);

    // *****************
    // init subcommand
    // *****************
    CLI::App* init = app.add_subcommand("init", "Initialize a set of subsequences to be used for embedding");

    int subseq_len;
    init->add_option("-l,--length", subseq_len, "Length (number of tokens) of the generated subsequences")
	->required();

    int token_len;
    init->add_option("-t,--token", token_len, "Length (number of characters) of a token")
	->required();

    int num_subseqs;
    init->add_option("-n,--number", num_subseqs, "Number of subsequences to generate, if input file(s) is provided, this is the number of subsequences to sample from *EACH* input sequence")
	->required();

    std::string alphabet_file;
    init->add_option("-a,--alphabet", alphabet_file, "File containing all permissible characters to be used")
	->default_val("../alphabets/DNA");
    
    std::vector<std::string> input_files;
    init->add_option("-i,--input", input_files, "Fasta file(s) to randomly sample subsequences from")
	->check(CLI::ExistingFile);

    std::string subseq_file;
    init->add_option("-o,--output", subseq_file, "File for storing the generated subsequences")
	->default_val("subsequences.txt");

    
    // *****************
    // sketch subcommand
    // *****************   
    CLI::App* sketch = app.add_subcommand("sketch", "Compute embedding of strings in the input file with the given subsequences");

    sketch->add_option("-s,--subsequences", subseq_file, "File containing the subsequences to be used for embedding")
	->required()
	->check(CLI::ExistingFile);

    sketch->add_option("-i,--input,fasta_files", input_files, "Fasta file(s) containing sequences to sketch")
	->required()
	->check(CLI::ExistingFile);

    
    // *****************
    // dist subcommand
    // *****************  
    CLI::App* dist = app.add_subcommand("dist", "Compute pairwise embedding distances between two embedding files");

    std::string embed_file1;
    dist->add_option("-a,--input1,embed_file1", embed_file1, "First file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string embed_file2;
    dist->add_option("-b,--input2,embed_file2", embed_file2, "Second file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string dist_file;
    dist->add_option("-o,--output", dist_file, "File for storing the embedding distances")
	->default_val("dist.rssebd-dist");

    
    // *****************
    // info subcommand
    // *****************   
    CLI::App* info = app.add_subcommand("info", "Show content of a binary embedding file");

    std::string embed_file;
    info->add_option("-i,--input,embed_file", embed_file, "Input embedding file")
	->required()
	->check(CLI::ExistingFile);

    
    // *****************
    // show subcommand
    // *****************   
    CLI::App* show = app.add_subcommand("show", "Show an embedding distance matrix stored in a binary file");

    show->add_option("-i,--input,dist_file", dist_file, "Input distance matrix file")
	->required()
	->check(CLI::ExistingFile);

    bool dist_to_stdout = true;
    show->add_flag("-o,--to-stdout,!-p,!--to-npy", dist_to_stdout, "Output the distance matrix to stdout (if -o) or to a npy file (if -p)");
    

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand(init))
    {
	gen_random_subsequences(subseq_len, token_len, num_subseqs,
				alphabet_file, input_files, subseq_file);
    }
    else if(app.got_subcommand(sketch))
    {
	compute_embeddings(subseq_file, input_files);
    }
    else if(app.got_subcommand(dist))
    {
	compute_distances(embed_file1, embed_file2, dist_file);
    }
    else if(app.got_subcommand(info))
    {
	show_embeddings(embed_file);
    }
    else if(app.got_subcommand(show))
    {
	show_distances(dist_file, dist_to_stdout);
    }
    
    return 0;
}

std::string change_file_ext(const std::string& file, const std::string& new_ext)
{
    std::size_t dot_pos = file.find_last_of('.');
    if(dot_pos == std::string::npos)
    {
	return file + '.' + new_ext;
    }
    else
    {
	return file.substr(0, dot_pos) + '.' + new_ext;
    }
}

void gen_random_subsequences(const int subseq_len,
			     const int token_len,
			     const int num_subseqs,
			     const std::string& alphabet_file,
			     const std::vector<std::string>& input_files,
			     const std::string& subseq_file)
{
    std::cout << "Generate random subsequences" << std::endl;
    std::cout << "subseq_len: " << subseq_len << std::endl;
    std::cout << "token_len: " << token_len << std::endl;
    std::cout << "num_subseqs: " << num_subseqs << std::endl;
    std::cout << "alphabet_file: " << alphabet_file << std::endl;
    std::cout << "input_files:";
    for(const std::string& s : input_files)
    {
	std::cout << " " << s;
    }
    std::cout << std::endl << "subseq_file: " << subseq_file
	      << std::endl << std::endl;

    std::cout << "Generating subsequences..." << std::endl;
    subsequences seqs(subseq_len, token_len);
    if(input_files.size() > 0)
    {
	seqs.gen_subsequences(input_files, num_subseqs);
    }
    else
    {
	std::ifstream fin(alphabet_file);
	if(!fin)
	{
	    // throw std::runtime_error("Could not open the file: " + alphabet_file);
	    std::cerr << "Error: could not open the file: "
		      << alphabet_file << std::endl;
	    std::exit(1);
	}
	std::string alphabet;
	if(!std::getline(fin, alphabet) || alphabet.empty())
	{
	    // throw std::runtime_error("Could not read alphabet from: " + alphabet_file);
	    std::cerr << "Error: could not read alphabet from: "
		      << alphabet_file << std::endl;
	    std::exit(1);
	}
	else
	{
	    std::cout << "Using alphabet: " << alphabet << std::endl;
	}

	fin.close();
	
	seqs.gen_subsequences(alphabet, num_subseqs);
    }

    seqs.save_subsequences(subseq_file);
    std::cout << "Generated " << seqs.size() << " subsequences, saved in "
	      << subseq_file << std::endl;
}

// Return the maximum number of consecutive tokens (starting from the
// leftmost one) in test that form a subsequence (of tokens) of seq
// using linear search.
int64_t longest_subsequence(const std::string& seq, const std::string& test,
			    int token_len)
{
    int result = 0;
    int64_t p = -1;

    for(int i = 0; i < test.size(); i += token_len)
    {
	p = seq.find(test.substr(i, token_len), p + 1);
	if(p == std::string::npos)
	{
	    break;
	}
	else
	{
	    result += 1;
	}
    }

    return result;
}

void compute_embeddings(const std::string& subseq_file,
			const std::vector<std::string>& input_files)
{
    std::cout << "Sketching" << std::endl << "input_files:";
    for(const std::string& s : input_files)
    {
	std::cout << " " << s;
    }
    std::cout << std::endl << "subseq_file: " << subseq_file
	      << std::endl << std::endl;

    subsequences subs(subseq_file);
    std::cout << "Loaded " << subs.size() << " subsequence(s), num_tokens: "
	      << subs.num_tokens << " token_len: "
	      << subs.token_len << std::endl;

    int num_subs = subs.size();
    std::string ext_name = "D" + std::to_string(num_subs) +
	".l" + std::to_string(subs.num_tokens) +
	".t" + std::to_string(subs.token_len) +
	".rssebd";
    
    for(const std::string& file : input_files)
    {	
	fasta_reader fin(file);
		
	std::vector<std::string> seqs;
	fin.read_all(seqs);
	size_t ct = seqs.size();

	std::cout << "Embedding " << ct << " sequence(s) in file: " << file << std::endl;

	Eigen::MatrixXi embeds(ct, num_subs);

#pragma omp parallel for default(shared)
	for(size_t i = 0; i < ct; ++i)
	{
	    for(int j = 0; j < num_subs; ++j)
	    {
		embeds(i, j) = longest_subsequence(seqs[i], subs.seqs[j], subs.token_len);
	    }
	}

	std::string out_file = change_file_ext(file, ext_name);
	rssebd_array::write_all(embeds, ct, num_subs, subs.num_tokens, out_file);
	
	std::cout << "Finished " << ct << " sequence(s), embedding wrote to file "
		  << out_file << std::endl;
    }
}


void compute_distances(const std::string& embed_file1,
		       const std::string& embed_file2,
		       const std::string& dist_file)
{
    std::cout << "embed_file1: " << embed_file1 << std::endl;
    std::cout << "embed_file2: " << embed_file2 << std::endl;
    std::cout << "dist_file: " << dist_file << std::endl << std::endl;

    std::cout << "Loading embeddings from the file: " << embed_file1 << std::endl;
    size_t num_embeds1;
    int embed_dim1;
    int num_tokens1;
    Eigen::MatrixXd embeds1 = rssebd_array::load_all(num_embeds1, embed_dim1, num_tokens1, true, false, embed_file1);
    std::cout << "Loaded " << num_embeds1 << " embeddings from "
	      << embed_file1 << ", dimension: " << embed_dim1 << std::endl;

    std::cout << "Loading embeddings from the file: " << embed_file2 << std::endl;
    size_t num_embeds2;
    int embed_dim2;
    int num_tokens2;
    Eigen::MatrixXd embeds2_tran = rssebd_array::load_all(num_embeds2, embed_dim2, num_tokens2, true, true, embed_file2);
    std::cout << "Loaded " << num_embeds2 << " embeddings from "
	      << embed_file2 << ", dimension: " << embed_dim2 << std::endl;

    if(embed_dim1 != embed_dim2)
    {
	std::cerr << "Error: embedding dimensions do not match." << std::endl;
	std::exit(1);
    }

    if(num_tokens1 != num_tokens2)
    {
	std::cerr << "Warning: max possible values in the embeddings are not consistent, #1: "
		  << num_tokens1 << ", #2: " << num_tokens2
		  << ". The results may not be meaningful." << std::endl;
    }

    std::cout << "Computing pairwise embedding distances..." << std::endl;
    // rssebd_array::pairwise_cos_dist(embeds1, embeds2, embed_dim1, dist_file);
    rssebd_array::pairwise_cos_dist(embeds1, embeds2_tran, dist_file);
    std::cout << num_embeds1 << "x" << num_embeds2
	      << " embedding distance matrix wrote to file: "
	      << dist_file << std::endl;

    // rssebd_array::free(embeds1);
    // rssebd_array::free(embeds2);
}


void show_embeddings(const std::string& embed_file)
{
    size_t num_embeds;
    int embed_dim;
    int num_tokens;

    std::cout << "Loading embeddings from the file: " << embed_file << std::endl;
    // rssebd_array::load(embeds, embed_dim, num_tokens, embed_file);
    Eigen::MatrixXd embeds = rssebd_array::load_all(num_embeds, embed_dim, num_tokens, false, false, embed_file);

    std::cout << "Embedding dimension: " << embed_dim << std::endl;
    std::cout << "Max possible value: " << num_tokens << std::endl;
    std::cout << "Number of embeddings: " << num_embeds << std::endl;   

    Eigen::IOFormat fmt(0, 0, " ", "\n");
    std::cout << embeds.format(fmt) << std::endl;
}


void show_distances(const std::string& dist_file, bool to_stdout)
{
    std::cout << "Loading distances from the file: " << dist_file << std::endl;
    rssebd_array::load_dist_matrix(dist_file, to_stdout);
}
