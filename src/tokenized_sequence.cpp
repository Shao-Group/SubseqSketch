#include "tokenized_sequence.hpp"

// TODO: consider using rolling hash if this is too slow
tokenized_sequence::tokenized_sequence(const std::string& seq, int token_len)
    : token_len(token_len)
{
    int64_t len = seq.size() - token_len + 1;
    for(int64_t i = 0; i < len; ++i)
    {
	auto result = index.emplace(seq.substr(i, token_len),
				    std::vector<int64_t>(1, i));
	// string exists in index, append index i
	if(!result.second)
	{
	    result.first->second.push_back(i);
	}
	
    }
}

int tokenized_sequence::longest_subsequence(const std::string& test) const
{
    int result = 0;
    int64_t p = -1;

    // assumes test.size is a multiple of token_len
    for(int i = 0; i < test.size(); i += token_len)
    {
	p = find(test.substr(i, token_len), p);
	if(p < 0)
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

int64_t tokenized_sequence::find(const std::string& token, int64_t st_pos) const
{
    auto it = index.find(token);
    if(it == index.end()) return -1;

    const std::vector<int64_t>& positions = it->second;

    int64_t result;
    
    if(positions.front() > st_pos) result =  positions[0];
    else if(positions.back() <= st_pos) result = -1;
    else // do binary search on positions to find smallest value > st_pos
    {
	int i = 1;
	int j = positions.size() - 1;
	while(i < j)
	{
	    int m = (i + j) >> 1;
	    if(positions[m] > st_pos)
	    {
		j = m;
	    }
	    else
	    {
		i = m + 1;
	    }
	}
	result = positions[i];
    }

    return result;
}
