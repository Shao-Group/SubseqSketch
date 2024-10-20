/*
  Part of SubseqEmbed.
  A naive fasta reader.
  By Ke @ Penn State
*/

#include "fasta_reader.hpp"
#include <sstream>

fasta_reader::fasta_reader(const std::string& file)
    : fin(file)
{
    if(!fin)
    {
	throw std::runtime_error("Could not open the file: " + file);
    }

    int header = fin.peek();
    if(header != '>')
    {
	throw std::runtime_error(file + " does not appear to be a valid fasta file");
    }
}

fasta_reader::~fasta_reader()
{
    fin.close();
}

bool fasta_reader::eof()
{
    return fin.eof();
}

std::string fasta_reader::next()
{
    assert(!eof());
#ifndef NDEBUG
    char header = fin.peek();
    assert(header == '>');
#endif

    // ignore the assumed header line
    fin.ignore(MAX_SIZE, '\n');
    std::string line;
    std::ostringstream oss;
    
    while(fin.peek() != '>' && std::getline(fin, line))
    {
	oss << line;
    }

    return oss.str();
}
