/*
  Part of SubseqEmbed.
  Representing a list of subsequences to be used for embedding.
  By Ke @ Penn State
*/

#include "subsequences.hpp"
#include "fasta_reader.hpp"
#include <fstream>
#include <sstream>
#include <random>

subsequences::subsequences(int subseq_len, int token_len)
    :token_len(token_len), num_tokens(subseq_len)
{}

subsequences::subsequences(const std::string& subseq_file)
{
    load_subsequences(subseq_file);
}

void subsequences::gen_subsequences(const std::string& alphabet, int num)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> dist(0, alphabet.size() - 1);

    int seq_len = token_len * num_tokens;

    seqs.reserve(num);
    
    std::string s;
    for(int i = 0; i < num; ++i)
    {
	s.reserve(seq_len);
	for(int j = 0; j < seq_len; ++j)
	{
	    s += alphabet[dist(generator)];
	}
	seqs.push_back(std::move(s));
    }
}

void subsequences::gen_subsequences(const std::vector<std::string>& input_files, int num_each)
{
    for(const std::string& file : input_files)
    {
	fasta_reader fin(file);
	while(!fin.eof())
	{
	    sample_subsequences(fin.next(), num_each);
	}
    }
}

void subsequences::save_subsequences(const std::string& subseq_file)
{
    std::ofstream fout(subseq_file);
    if(!fout)
    {
	throw std::runtime_error("Could not write to file: " + subseq_file);
    }

    fout << seqs.size() << " " << num_tokens << " " << token_len << std::endl;
    for(const std::string& s : seqs)
    {
	fout << s << std::endl;
    }
    fout.close();
}

void subsequences::load_subsequences(const std::string& subseq_file)
{
    std::ifstream fin(subseq_file);
    if(!fin)
    {
	throw std::runtime_error("Could not open file: " + subseq_file);
    }

    int num_seqs;
    fin >> num_seqs >> num_tokens >> token_len;
    fin.ignore(MAX_SIZE, '\n');

    seqs.reserve(num_seqs);
    std::string s;
    while(std::getline(fin, s))
    {
	seqs.push_back(std::move(s));
    }

    fin.close();
}

void subsequences::sample_subsequences(const std::string& reference, int num)
{
    std::random_device rd;
    std::mt19937 generator(rd());

    size_t part_len = reference.size() / num_tokens;
    assert(part_len > token_len);
    
    std::uniform_int_distribution<> dist(0, part_len - token_len);

    seqs.reserve(seqs.size() + num);

    std::string s;
    int seq_len = token_len * num_tokens;
    
    for(int i = 0; i < num; ++i)
    {
	s.reserve(seq_len);
	size_t p = 0;
	for(int j = 0; j < num_tokens; ++j)
	{
	    s.append(reference, p + dist(generator), token_len);
	    p += part_len;
	}
	seqs.push_back(std::move(s));
    }
}

std::size_t subsequences::size() const
{
    return seqs.size();
}
