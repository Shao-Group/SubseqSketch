/*
  Part of SubseqSketch.
  Representing a list of subsequences to be used for sketching.
  By Ke @ Penn State
*/

#ifndef __SUBSEQUENCES_H__
#define __SUBSEQUENCES_H__

#include <vector>
#include <string>

class subsequences
{
public:
    int token_len;
    int num_tokens;
    std::vector<std::string> seqs;

    // Initialize an empty list of subsequences.
    subsequences(int subseq_len, int token_len);
    // Load a list of subsequences from the given file.
    subsequences(const std::string& subseq_file);

    // Generate num random subsequences on the given alphabet.
    void gen_subsequences(const std::string& alphabet, int num);
    // Sample num_each random subsequences from each of the sequence in each of
    // the given fasta files.
    void gen_subsequences(const std::vector<std::string>& input_files, int num_each);

    // Write the subsequences to a file.
    // First line contains three numbers for number of subsequences,
    // number of tokens per subsequence, and length of a token.
    // Following lines each contains a subsequence.
    void save_subsequences(const std::string& subseq_file);

    std::size_t size() const;

private:
    void load_subsequences(const std::string& subseq_file);

    // Sample num subsequences from the given reference.
    // The reference is first split into num_tokens parts, one random token
    // is then sampled from each part to form a sampled subsequence.
    void sample_subsequences(const std::string& reference, int num);
};


#endif
