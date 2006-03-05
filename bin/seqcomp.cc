#include <iostream>
#include <unistd.h>

#include "paircomp.hh"

#define BIN_VERSION "1.2"

// function for reading in FASTA sequences
std::string read_fasta_seq(char * filename);

// length cutoff for input sequences
#define LENGTH_CUTOFF 300000

//
// main function
//

int main(int argc, char * argv[])
{
  char * file_in1 = NULL, * file_in2 = NULL, * file_out = NULL;
  unsigned int windowsize = 0;
  unsigned int ithreshold = 0;
  float threshold = 1.0;
  bool do_exit = false;

  // Parse options.
  while (1) {
    char ch = getopt(argc, argv, "v");
    if (ch == 'v') {
      std::cout << "paircomp toolkit version " << VERSION \
		<< "; seqcomp compatibility binary v" << BIN_VERSION << "\n";
      do_exit = true;
    }
    else if (ch == -1) {		// end of options
      break;
    }
  }

  if (do_exit) {
    exit(0);
  }

  if (optind != argc - 6) {
    std::cerr << "Usage:\n\t" << argv[0] \
	      << " seq1 seq2 windowsize ithreshold 'noxml' out\n\n";
    exit(-1);
  }

  file_in1 = argv[optind];
  file_in2 = argv[optind + 1];
  windowsize = atoi(argv[optind + 2]);
  ithreshold = atoi(argv[optind + 3]);
  file_out = argv[optind + 5];

  // Verify that windowsize & threshold are reasonable.
  if (windowsize < 10 || windowsize > 500) {
    std::cerr << "ERROR: windowsize must be 10 <= w < 500.\n";
    exit(-1);
  }

  if (ithreshold > windowsize) {
    std::cerr << "ERROR: ithreshold must be <= the windowsize.\n";
    exit(-1);
  }

  threshold = (float) ithreshold / (float) windowsize;

  // Check that we can open the output file for reading.
  FILE * fp = fopen(file_out, "w");
  if (fp == NULL ){
    std::cerr << "ERROR: cannot open file '" << file_out << "' for writing.\n";
    exit(-1);
  } else {
    fclose(fp);
  }

  // Read in the sequences.
  std::string seq1 = read_fasta_seq(file_in1);
  std::string seq2 = read_fasta_seq(file_in2);

  if (seq1.length() > LENGTH_CUTOFF || seq2.length() > LENGTH_CUTOFF) {
    std::cerr << "ERROR: sequence too long; length cutoff is "
	      << LENGTH_CUTOFF << "\n";
    exit(-1);
  }

  std::cout << "Doing a paircomp with windowsize " << windowsize \
	    << " and threshold " << threshold << ".\n";

  // Do paircomp.
  paircomp::ImmutableComparison * cmp;
  cmp = paircomp::rolling_nxn_comparison(seq1, seq2, windowsize, threshold);

  // Write to output file.
  cmp->save_as_seqcomp(file_out);

  delete cmp;

  exit(0);
}

// Read in a FASTA sequence from a file.

std::string read_fasta_seq(char * filename)
{
  FILE * fp = fopen(filename, "r");
  if (fp == NULL) {
    std::cerr << "ERROR: Cannot open filename '" << filename \
	      << "' for reading.\n";
    exit(-1);
  }

  int filesize = 0;
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  rewind(fp);

  if (filesize <= 0) {
    std::cerr << "ERROR: no data in file '" << filename << "'???\n";
    exit(-1);
  }

  // allocate a new buffer of that size.
  char * buf = new char[filesize+1];

  // Read in the file.
  long items_read = fread(buf, filesize, 1, fp);
  if (items_read != 1) {
    std::cerr << "ERROR: short read on file '" << filename << "'.\n";
    exit(-1);
  }
  buf[filesize] = 0;

  std::string contents(buf);

  delete buf;

  if (contents[0] != '>') {
    std::cerr << "ERROR: contents of file '" << filename \
	      << "' does not begin with a '>'.\n";
    exit(-1);
  }

  // Separate out the sequence.
  std::string::size_type pos = contents.find("\n");
  
  std::string name = std::string(contents, 1, pos - 1);
  std::string seq = std::string(contents, pos + 1, contents.length());

  // Iterate through & remove whitespace.
  pos = seq.find("\n");
  while (pos != std::string::npos) {
    seq.erase(pos, 1);
    pos = seq.find("\n");
  }

  // Make sure there's only one sequence.
  if (seq.find(">") != std::string::npos) {
    std::cerr << "ERROR: file '" << filename \
	      << "' contains more than one sequence.\n";
    exit(-1);
  }

  std::cout << "Loaded sequence '" << name << "' from file '" \
	    << filename << "'; length " << seq.length() << ".\n";

  return seq;
}
