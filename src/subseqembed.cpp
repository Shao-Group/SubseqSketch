/*
  Edit distance embedding by random subsequences.
  By Ke @ Penn State
*/

#include <iostream>
#include <vector>
#include "CLI11.hpp"

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
    init->add_option("-n,--number", num_subseqs, "Number of subsequences to generate, if input file(s) is provided, this is the number of subsequences to sample from EACH input sequence")
	->required();

    std::string alphabet_file;
    init->add_option("-a,--alphabet", alphabet_file, "File containing all permissible characters to be used")
	->default_val("alphabets/DNA");
    
    std::vector<std::string> input_files;
    init->add_option("-i,--input", input_files, "Fasta files to randomly sample subsequences from")
	->each(CLI::ExistingFile);

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

    sketch->add_option("-i,--input", input_files, "Fasta files containing sequences to sketch")
	->required()
	->each(CLI::ExistingFile);


    // *****************
    // dist subcommand
    // *****************  
    CLI::App* dist = app.add_subcommand("dist", "Compute pairwise embedding distances between two embedding files");

    std::string embed_file1;
    dist->add_option("-a,--input1", embed_file1, "First file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string embed_file2;
    dist->add_option("-b,--input2", embed_file2, "Second file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string dist_file;
    dist->add_option("-o,--output", dist_file, "File for storing the embedding distances")
	->required();

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
    
    return 0;
}

void gen_random_subsequences(const int subseq_len,
			     const int token_len,
			     const int num_subseqs,
			     const std::string& alphabet_file,
			     const std::vector<std::string>& input_files,
			     const std::string& subseq_file)
{
    std::cout << "subseq_len: " << subseq_len << std::endl;
    std::cout << "token_len: " << token_len << std::endl;
    std::cout << "num_subseqs: " << num_subseqs << std::endl;
    std::cout << "alphabet_file: " << alphabet_file << std::endl;
    std::cout << "input_files:";
    for(const std::string& s : input_files)
    {
	std::cout << " " << s;
    }
    std::cout << std::endl << "subseq_file: " << subseq_file << std::endl;
}

void compute_embeddings(const std::string& subseq_file,
			const std::vector<std::string>& input_files)
{
    std::cout << "input_files:";
    for(const std::string& s : input_files)
    {
	std::cout << " " << s;
    }
    std::cout << std::endl << "subseq_file: " << subseq_file << std::endl;
}

void compute_distances(const std::string& embed_file1,
		       const std::string& embed_file2,
		       const std::string& dist_file)
{
    std::cout << "embed_file1: " << embed_file1 << std::endl;
    std::cout << "embed_file2: " << embed_file2 << std::endl;
    std::cout << "dist_file: " << dist_file << std::endl;
}
