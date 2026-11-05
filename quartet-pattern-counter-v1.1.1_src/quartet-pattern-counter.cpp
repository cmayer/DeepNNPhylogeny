// Copyright (c) 2022, Christoph Mayer and Nikita Kulikov

// Redistribution and use in source and binary forms, 
// with or without modification, are permitted provided
// that the following conditions are met:

// Redistributions of source code must retain the above
// copyright notice, this list of conditions and the
// following disclaimer.

// Redistributions in binary form must reproduce the
// above copyright notice, this list of conditions and
// the following disclaimer in the documentation and/or
// other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
// CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>
#include "symbol-cartesian-product_string.h"
#include "numpy.h"

using namespace std;

void read_fasta_file_into_map(const char *filename, map<string, string> &fastamap, vector<string> &namevec)
{
  ifstream is;
  is.open(filename);
  string line, lastseqstring, lastname;

  fastamap.clear();

  if (is.fail())
  {
    return;
  }

  while (getline(is,line))
  {
    if (line[0] == '>')
    {
      if (!lastseqstring.empty())
      {
	fastamap[lastname] = lastseqstring;
	namevec.push_back(lastname);
      }
      lastname = line.substr(1, string::npos);
      lastseqstring.clear();
    }
    else
    {
      line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
      lastseqstring.append(line);
    }
  }
  if (!lastseqstring.empty())
  {
    fastamap[lastname] = lastseqstring;
    namevec.push_back(lastname);
  }
  is.close();
}


// Returns the number of sites that were counted. Sites with a symbol that is not part of the
// alphabet (gaps, N, ambiguity codes, ...) in any of the four sequences are skipped.
// Lower case symbols are counted as their upper case equivalent.
unsigned get_pattern_map(const char **four_sequences,
		     unsigned N, map<string, unsigned> &pattern_counter, bool dataType_protein)
{
  unsigned num_valid_sites = 0;
  vector<string> all_patterns;
  cartesian_product  cp;
  if (dataType_protein)
    cp.operator()(4, "ARNDCQEGHILKMFPSTWYV", all_patterns);
  else
    cp.operator()(4, "ACGT", all_patterns);

  string pat = "    ";

  for (unsigned i=0; i<all_patterns.size(); ++i)
  {
    pattern_counter[all_patterns[i]] = 0;
  }
  
  for (unsigned pos=0; pos<N; ++pos)
  {
    pat[0] = toupper((unsigned char)four_sequences[0][pos]);
    pat[1] = toupper((unsigned char)four_sequences[1][pos]);
    pat[2] = toupper((unsigned char)four_sequences[2][pos]);
    pat[3] = toupper((unsigned char)four_sequences[3][pos]);

    // Only count patterns of valid symbols. Inserting other patterns into the map
    // would shift the output vector, which is written in map order.
    map<string, unsigned>::iterator it = pattern_counter.find(pat);
    if (it == pattern_counter.end())
      continue;
    ++(it->second);
    ++num_valid_sites;
  }
  return num_valid_sites;
}


void write_pattern_fequencies_to_file(map<string, unsigned> &pattern_counter,
				      unsigned num_sites, const char *outfilename,
				      bool dataType_protein)
{
  map<string, unsigned>::iterator it_pat, it_pat_end;
  it_pat     = pattern_counter.begin();
  it_pat_end = pattern_counter.end();
  vector<float> pattern_frequencies;

  if (dataType_protein)
    pattern_frequencies.reserve(160000);
  else
    pattern_frequencies.reserve(256);

  //  cerr << "pattern_counter size: " << pattern_counter.size() << endl;

  while (it_pat != it_pat_end)
  {
    pattern_frequencies.push_back((double)it_pat->second/(double)num_sites);
    //    cout << it_pat->first << " " << (double)it_pat->second/(double)num_sites << endl;
    ++it_pat;
  }
  if (dataType_protein)
    aoba::SaveArrayAsNumpy(outfilename, 160000, &pattern_frequencies[0]);
  else
    aoba::SaveArrayAsNumpy(outfilename, 256, &pattern_frequencies[0]);
}


int main(int argc, char **argv)
{
  bool dataType_protein = false;
  char *filename;
  char *outfilename;
  const char *four_sequences[4];

  if (argc < 3 || argc > 4)
  {
    cerr << "quartet-pattern-counter v1.1.1\nUsage: quartet-pattern-counter-v1.1 [-p] infile.fas outfile.npy\n";
    exit(-1);
  }

  int i,j;
  for (i=1, j=1; i<argc; ++i)
  {
    if (string(argv[i]) == string("-p") )
    {
      dataType_protein = true;
    }
    else if (j==1)
    {
      filename = argv[i];
      ++j;
    }
    else if (j == 2)
    {
      outfilename = argv[i];
    }
  }

  map<string, string> fasta_map_string;
  map<string, unsigned> pattern_counter;
  vector<string> seq_names;
  unsigned len = 0;
    
  read_fasta_file_into_map(filename, fasta_map_string, seq_names);

  if (seq_names.size() != 4 || fasta_map_string.size() != 4)
  {
    cerr << "Error: The input file " << filename << " must contain exactly four sequences with distinct names. Found: "
	 << seq_names.size() << " sequences.\n";
    exit(-1);
  }

  for (unsigned i=0; i<seq_names.size(); ++i)
  {
      //       cout << ">>" << seq_names[i] << '\n'
      //            << fasta_map_string[seq_names[i]] << '\n';
    string &tmp = fasta_map_string[seq_names[i]];
    if (i==0)
      len = tmp.length();
    else if (tmp.length() != len)
    {
      cerr << "Error: The sequences in the input file " << filename << " do not have equal lengths.\n";
      exit(-1);
    }
    four_sequences[i] = tmp.c_str();
  }
  unsigned num_valid_sites = get_pattern_map(four_sequences, len, pattern_counter, dataType_protein);
  if (num_valid_sites == 0)
  {
    cerr << "Error: The input file " << filename << " contains no site without gaps or ambiguous symbols.\n";
    exit(-1);
  }
  write_pattern_fequencies_to_file(pattern_counter, num_valid_sites, outfilename, dataType_protein);

  exit(0);


}
