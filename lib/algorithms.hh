
#ifndef ALGORITHMS_HH
#define ALGORITHMS_HH

#include <string>
#include "MutableComparison.hh"
#include "ImmutableComparison.hh"

namespace paircomp {
  ImmutableComparison * simple_nxn_comparison(const std::string seq1,
					      const std::string seq2,
					      unsigned int windowsize,
					      float threshold);

  ImmutableComparison * rolling_nxn_comparison(const std::string seq1,
					       const std::string seq2,
					       unsigned int windowsize,
					       float threshold);

  //
  // Some useful utility functions.
  //

  // transform sequence to uppercase.
  void _preprocess(std::string &seq);
  std::string _reverse_complement(const std::string seq);

  bool _match_window(const char * top, unsigned int start_top,
		     const char * bot, unsigned int start_bot,
		     const unsigned int windowsize,
		     const unsigned int threshold,
		     unsigned int &matches);

  void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const int orientation);

  void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const char * bot_seq_rev);

  unsigned int paircomp::_n_matching(const char * a, unsigned int a_start,
				     const char * b, unsigned int b_start,
				     unsigned int w);

  void _check_parameters(std::string top, std::string top,
			 unsigned int, float);
}

#endif // ALGORITHMS_HH
