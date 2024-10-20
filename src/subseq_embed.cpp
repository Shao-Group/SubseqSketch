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
#include "tokenized_sequence.hpp"
#include "rssebd_array.hpp"
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

void show_embeddings(const std::string& embed_file);

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
    dist->add_option("-a,--input1", embed_file1, "First file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string embed_file2;
    dist->add_option("-b,--input2", embed_file2, "Second file of embeddings")
	->required()
	->check(CLI::ExistingFile);

    std::string dist_file;
    dist->add_option("-o,--output", dist_file, "File for storing the embedding distances")
	->default_val("dist.txt");

    
    // *****************
    // info subcommand
    // *****************   
    CLI::App* info = app.add_subcommand("info", "Show content of a binary embedding file");

    std::string embed_file;
    info->add_option("-i,--input,embed_file", embed_file, "Input embedding file")
	->required()
	->check(CLI::ExistingFile);


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
    std::cout << std::endl << "subseq_file: " << subseq_file << std::endl;

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
	    throw std::runtime_error("Could not open the file: " + alphabet_file);
	}
	std::string alphabet;
	if(!std::getline(fin, alphabet) || alphabet.empty())
	{
	    throw std::runtime_error("Could not read alphabet from: " + alphabet_file);
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


void compute_embeddings(const std::string& subseq_file,
			const std::vector<std::string>& input_files)
{
    std::cout << "Sketching" << std::endl << "input_files:";
    for(const std::string& s : input_files)
    {
	std::cout << " " << s;
    }
    std::cout << std::endl << "subseq_file: " << subseq_file << std::endl;

    subsequences subs(subseq_file);
    std::cout << "Loaded " << subs.size() << " subsequence(s), num_tokens: "
	      << subs.num_tokens << " token_len: "
	      << subs.token_len << std::endl;

    int num_subs = subs.size();
    int* embed = new int[num_subs];
    
    for(const std::string& file : input_files)
    {
	std::cout << "Embedding sequence(s) in file: " << file << std::endl;
	
	fasta_reader fin(file);
	std::string out_file = change_file_ext(file, "rssebd");
	std::ofstream fout(out_file, std::ios::binary);
	if(!fout)
	{
	    throw std::runtime_error("Could not write to file: " + out_file);
	}
	
	int ct = 0;
	while(!fin.eof())
	{
	    tokenized_sequence s_index(fin.next(), subs.token_len);
	    
	    for(int i = 0; i < num_subs; ++i)
	    {
		embed[i] = s_index.longest_subsequence(subs.seqs[i]);
	    }

	    rssebd_array::write(embed, num_subs, subs.num_tokens, fout);
	    ++ct;
	}

	std::cout << "Finished " << ct << " sequence(s), embedding wrote to file "
		  << out_file << std::endl;
    }

    delete[] embed;
}


void compute_distances(const std::string& embed_file1,
		       const std::string& embed_file2,
		       const std::string& dist_file)
{
    std::cout << "embed_file1: " << embed_file1 << std::endl;
    std::cout << "embed_file2: " << embed_file2 << std::endl;
    std::cout << "dist_file: " << dist_file << std::endl;
}


void show_embeddings(const std::string& embed_file)
{
    std::vector<int*> embeds;
    int embed_dim;
    int num_tokens;

    rssebd_array::load(embeds, embed_dim, num_tokens, embed_file);

    std::cout << "Embedding dimension: " << embed_dim << std::endl;
    std::cout << "Max possible value: " << num_tokens << std::endl;
    std::cout << "Number of embeddings: " << embeds.size() << std::endl;   
    
    for(int* cur : embeds)
    {
	std::cout << cur[0];
	for(int i = 1; i < embed_dim; ++i)
	{
	    std::cout << " " << cur[i];
	}
	std::cout << std::endl;
	delete[] cur;
    }
}
