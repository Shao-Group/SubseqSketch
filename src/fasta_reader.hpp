/*
  Part of SubseqEmbed.
  A naive fasta reader.
  By Ke @ Penn State
*/

#ifndef __FASTA_READER_H__
#define __FASTA_READER_H__

#include <fstream>
#include <string>
#include <limits>
#include <iostream>
#include <cassert>
 
constexpr auto MAX_SIZE = std::numeric_limits<std::streamsize>::max();

class fasta_reader
{
public:
    fasta_reader(const std::string& file);
    ~fasta_reader();

    bool eof();
    // Read the next sequence in file. A sequence is assumed to be preceded by
    // a header line starts with '>', the header line is then ignored. Following
    // lines until the next header or eof are concatenated and returned.
    std::string next();
    
private:
    std::ifstream fin;
};

#endif
