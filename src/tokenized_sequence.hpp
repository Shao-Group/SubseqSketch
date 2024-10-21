/*
  Part of SubseqEmbed.
  Index of a string for fast (tokenized) subsequence searching.
  By Ke @ Penn State
*/

#ifndef __TOKENIZED_SEQUENCE_H__
#define __TOKENIZED_SEQUENCE_H__

#include <unordered_map>
#include <vector>
#include <string>
#include <cstdint>

class tokenized_sequence
{
public:
    tokenized_sequence(const std::string& seq, int token_len);
    // Return the maximum number of consecutive tokens (starting from the
    // leftmost one) in test that form a subsequence (of tokens) of this
    // underlying sequence.
    int longest_subsequence(const std::string& test) const;
private:
    int token_len;
    // For each token appears in the underlying sequence, the indices of its
    // occurrences are stored in ascending order in the corresponding vector.
    std::unordered_map<std::string, std::vector<int64_t>> index;

    // Search if token appears in the underlying string starting from st_pos+1.
    // If found, return the beginning index of that occurrence, otherwise
    // return -1.
    int64_t find(const std::string& token, int64_t st_pos) const;
};

#endif
