#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
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

/*
void read_fasta_file_into_map(const char *filename, map<faststring, faststring> &fastamap, vector<faststring> &namevec)
{
  ifstream is;
  is.open(filename);
  faststring line, lastseqstring, lastname;

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
      lastname = line;
      lastseqstring.clear();
      lastname = line.substr(1, faststring::npos);
    }
    else
      lastseqstring.append(line);
  }
  if (!lastseqstring.empty())
  {
    fastamap[lastname] = lastseqstring;
    namevec.push_back(lastname);
  }
}


void read_fasta_file_into_map(const char *filename,   CSequences2 **pseqs)
{
  CFile file;

  file.ffopen(filename);

  if (file.fail())
  {
    cerr << "Could not open specified file: " << filename << endl;
    exit(-1);
  }

  *pseqs = new CSequences2(CSequence_Mol::molecular);

  CSequence_Mol::processing_flag pflag;
  pflag  = CSequence_Mol::processing_flag(0);
  int error = (*pseqs)->read_from_Fasta_File(file, pflag, 0, -1, false);
}
*/

void get_pattern_map(const char **four_sequences,
		     unsigned N, map<string, unsigned> &pattern_counter, bool dataType_protein)
{
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
    pat[0] = four_sequences[0][pos];
    pat[1] = four_sequences[1][pos];
    pat[2] = four_sequences[2][pos];
    pat[3] = four_sequences[3][pos];
    
    //    cerr << "(1) Pattern increment: " << pat << " " << pattern_counter[pat] << endl; 
    ++pattern_counter[pat];
    //    cerr << "(2) Pattern increment: " << pat << " " << pattern_counter[pat] << endl; 
  }
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
    cerr << "Usage: pattern_counter_quartett_trees [-p]infile.fas outfile.npy\n";
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

  for (unsigned i=0; i<seq_names.size(); ++i)
  {
      //       cout << ">>" << seq_names[i] << '\n'
      //            << fasta_map_string[seq_names[i]] << '\n';
    string &tmp = fasta_map_string[seq_names[i]];
    if (i==0)
      len = tmp.length();
    four_sequences[i] = tmp.c_str();
  }
  get_pattern_map(four_sequences, len, pattern_counter, dataType_protein);
  write_pattern_fequencies_to_file(pattern_counter, len, outfilename, dataType_protein); 

  /*
  if (stringType == 1)
  {
    map<faststring, faststring> fasta_map_string;
    map<string, unsigned> pattern_counter;
    vector<faststring> seq_names;
    unsigned len = 0;
    
    read_fasta_file_into_map(filename, fasta_map_string, seq_names);

    for (unsigned i=0; i<seq_names.size(); ++i)
    {
       //       cout << ">>" << seq_names[i] << '\n' << fasta_map_string[seq_names[i]] << '\n';
      faststring &tmp = fasta_map_string[seq_names[i]];
      if (i==0)
	len = tmp.length();
      four_sequences[i] = tmp.c_str();
    }
  }

  if (stringType == 2)
  {
    map<faststring, faststring> fasta_map_string;
    map<string, unsigned> pattern_counter;
    vector<faststring> seq_names;
    
    CSequences2 *seqs;
    read_fasta_file_into_map(filename, &seqs);

    for (unsigned i=0; i<seqs->GetTaxaNum(); ++i)
    {
      //       cout << ">>" << seq_names[i] << '\n' << fasta_map_string[seq_names[i]] << '\n';
      seqs->get_Seq_Data(i);
    }
  }
  */
  exit(0);


}
