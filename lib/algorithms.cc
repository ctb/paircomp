// Function implementations for different 'seqcomp'-style algorithms.
//
// See README.txt for license & copyright information.
/*!
*   \file algorithms.cc
*/
#include <string>
#include <ctype.h>

#include "algorithms.hh"

using namespace paircomp;

/*! \def MATCH(top_s, top_i, bot_s, bot_i)
*   \brief Macro checks to see if characters in top_s at index top_i and bot_s at index bit_i match. 
*   \param top_s Character array for the top DNA/RNA/Protein sequence.
*   \param top_i Index for the top sequence.
*   \param bot_s Character array for the bottom DNA/RNA/Protein sequence.
*   \param bot_i Index for the bottom sequence.
*/
#define MATCH(top_s, top_i, bot_s, bot_i) ((top_s)[top_i] == (bot_s)[bot_i] &&\
	      (top_s)[top_i] != 'N' && \
	      (top_s)[top_i] != 'X')

#define TEST

bool paircomp::_match_window(const char * top,
			     const unsigned int start_top,
			     const char * bot,
			     const unsigned int start_bot,
			     const unsigned int w, const unsigned int th,
			     unsigned int &m)
{
  m = 0;

  if (start_top + w > strlen(top) || start_bot + w > strlen(bot)) {
    assert(0);
  }

  for (unsigned int k = 0; k < w; k++) {
    if (MATCH(top, start_top + k, bot, start_bot + k)) {
      m++;
    }
  }
  if (m >= th) {
    return true;
  }
  return false;
}

unsigned int paircomp::_n_matching(const char * a, unsigned int a_start,
				   const char * b, unsigned int b_start,
				   unsigned int w)
{
  unsigned int m = 0;
  for (unsigned int i = 0; i < w; i++) {
    if (MATCH(a, a_start + i, b, b_start +i)) {
      m++;
    }
  }
  return m;
}

//
// simple_nxn_comparison -- does a simple O(NxMxW) comparison, with
// three loops.  The code's pretty simple...
//

ImmutableComparison * paircomp::simple_nxn_comparison(const std::string seq1,
						      const std::string seq2,
						      unsigned int windowsize,
						      float threshold)
{
  // transform the strings appropriately: convert to uppercase, etc.
  std::string top_seq = seq1;
  std::string bot_seq = seq2;

  _preprocess(top_seq);
  _preprocess(bot_seq);

  _check_parameters(top_seq, bot_seq, windowsize, threshold);

  std::string seq2_rev = _reverse_complement(bot_seq);
  
  const char * top = top_seq.c_str();
  const char * bot = bot_seq.c_str();
  const char * bot_rev = seq2_rev.c_str();

  unsigned int min_matches = (unsigned int) (threshold * (float)windowsize + .5);

  const unsigned int top_end = top_seq.length() - windowsize + 1;
  const unsigned int top_len = top_seq.length();
  const unsigned int bot_end = bot_seq.length() - windowsize + 1;
  const unsigned int bot_len = bot_seq.length();

  // Create an empty mutable comparison.
  MutableComparison cmp(top_len, bot_len, windowsize);

  //
  // Do the brutal n^2 forward comparison.
  //

  for (unsigned int i = 0; i < top_end; i++) {
    for (unsigned int j = 0; j < bot_end; j++) {
      unsigned int m = 0;

      // should we record this window?
      if (_match_window(top, i, bot, j, windowsize, min_matches, m)) {
	Match * match = new Match(i, j, windowsize, m, 1);
	cmp.add_match(match);
      }
    }
  }

  //
  // Do the brutal n^2 reverse comparison.
  //

  for (unsigned int i = 0; i < top_end; i++) {
    for (unsigned int j = 0; j < bot_end; j++) {
      unsigned int m = 0;

      // should we record this window?
      if (_match_window(top, i, bot_rev, j, windowsize, min_matches, m)) {
	Match * match = new Match(i, j, windowsize, m, 1);
	match->reverse_bottom(bot_len);
	cmp.add_match(match);
      }
    }
  }

  // Return a more useful immutable comparison now that we're done.
  return new ImmutableComparison(&cmp);
}

////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////

//
// rolling_nxn_comparison(const std::string seq1, const std::string seq2,
//

ImmutableComparison * paircomp::rolling_nxn_comparison(const std::string seq1,
						       const std::string seq2,
						       unsigned int windowsize,
						       float threshold)
{
  // transform the strings appropriately: convert to uppercase, etc.
  std::string top_seq = seq1;
  std::string bot_seq = seq2;

  _preprocess(top_seq);
  _preprocess(bot_seq);

  _check_parameters(top_seq, bot_seq, windowsize, threshold);
  
  std::string seq2_rev = _reverse_complement(bot_seq);
  
  const char * top = top_seq.c_str();
  const char * bot = bot_seq.c_str();
  const char * bot_rev = seq2_rev.c_str();

  unsigned int min_matches = (unsigned int) (threshold * (float)windowsize + .5);

  const unsigned int top_len = top_seq.length();
  const unsigned int top_end = top_seq.length() - windowsize + 1;
  const unsigned int bot_len = bot_seq.length();
  const unsigned int bot_end = bot_seq.length() - windowsize + 1;

  // Create an empty mutable comparison.
  MutableComparison cmp(top_len, bot_len, windowsize);

  // Do the comparison...

  if (top_end == bot_end) { // sequences are same length
    // first do the zeroth offset.
    _compare_stretch(cmp, windowsize, min_matches,
		     top, top_end, 0,
		     bot, bot_end, 0,
		     bot_rev);

    for (unsigned int offset = 1; offset < top_end; offset++) {
      // do first part
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, 0,
		       bot, bot_end, offset,
		       bot_rev);

      // do second part
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, top_end - offset,
		       bot, bot_end, 0,
		       bot_rev);
    }
  } else if (top_end > bot_end) { // top sequence > bottom sequence
    for (unsigned int offset = 0; offset < top_end; offset++) {
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, offset,
		       bot, bot_end, 0,
		       bot_rev);
    }
    for (unsigned int offset = 1; offset < bot_end; offset++) {
      // do beginning of top
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, 0,
		       bot, bot_end, offset,
		       bot_rev);
    }
  } else if (top_end < bot_end) { // top sequence < bottom sequence
    for (unsigned int offset = 0; offset < bot_end; offset++) {
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, 0,
		       bot, bot_end, offset,
		       bot_rev);
    }
    for (unsigned int offset = 1; offset < top_end; offset++) {
      // do beginning of bot
      _compare_stretch(cmp, windowsize, min_matches,
		       top, top_end, offset,
		       bot, bot_end, 0,
		       bot_rev);
    }
    
  } else {
    assert(0);			// how can this happen?  ehh, leave it...
  }

  // Return a more useful immutable comparison now that we're done.
  return new ImmutableComparison(&cmp);
}

//
// Some useful utilities.
//

void paircomp::_preprocess(std::string &seq)
{
  // transform to uppercase
  transform(seq.begin(), seq.end(), seq.begin(), toupper);

  // make sure seq contains valid characters only
  std::string::size_type pos = seq.find_first_not_of("ACGTNX");

  if (pos != std::string::npos) {
    std::string msg = "error, invalid character '";
    msg += seq[pos];
    msg += "' (";
    msg += (int) seq[pos];
    msg += ") in sequence.";

    throw paircomp_exception(msg);
  }
}

static char _complement(char ch)
{
  switch(ch) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'N':
    return 'N';
  case 'X':
    return 'X';
  }
  throw paircomp_exception("unknown character");
  return 0;			// avoid compiler warning ;)
}

//
// Build reverse complement.
//

std::string paircomp::_reverse_complement(const std::string seq)
{
  std::string r(seq);
  reverse(r.begin(), r.end());	// first build reverse

  // then build complement
  transform(r.begin(), r.end(), r.begin(), _complement);

  return r;
}

//
// _compare_stretch
//

void paircomp::_compare_stretch(MutableComparison &cmp,
				unsigned int windowsize,
				unsigned int min_matches,
				const char * top_seq, unsigned int top_end,
				unsigned int top_offset,
				const char * bot_seq, unsigned int bot_end,
				unsigned int bot_offset,
				const char * bot_seq_rev)
{
  _compare_stretch(cmp, windowsize, min_matches,
		   top_seq, top_end, top_offset,
		   bot_seq, bot_end, bot_offset,
		   1);
  _compare_stretch(cmp, windowsize, min_matches,
		   top_seq, top_end, top_offset,
		   bot_seq_rev, bot_end, bot_offset,
		   -1);
}

void paircomp::_compare_stretch(MutableComparison &cmp,
				unsigned int windowsize,
				unsigned int min_matches,
				const char * top_seq, unsigned int top_end,
				unsigned int top_offset,
				const char * bot_seq, unsigned int bot_end,
				unsigned int bot_offset,
				const int orientation)
{
  assert(orientation == -1 || orientation == 1);
  Match * match;

  unsigned int i, k, matches;
  bool first_was_match;

  const unsigned int window_end = windowsize - 1;
  const unsigned int bot_len = strlen(bot_seq);

  assert (top_offset < top_end);
  assert (bot_offset < bot_end);
  assert (strlen(top_seq) == top_end + window_end);
  assert (strlen(bot_seq) == bot_end + window_end);

  // seed the search.
  matches = 0;
  for (k = 0; k < windowsize; k++) {
    if (MATCH(top_seq + top_offset, k, bot_seq, k + bot_offset)) {
      matches++;
    }
  }

  // should we record this window?
  if (matches >= min_matches) {
    match = new Match(top_offset, bot_offset, windowsize, matches, 1);
    if (orientation == -1) { match->reverse_bottom(bot_len); }

    cmp.add_match(match);
  }

  // was the first bp a match?
  first_was_match = MATCH(top_seq, top_offset, bot_seq, bot_offset);

  //
  // now iterate.
  //

  const unsigned int top_iterate_end = top_end - top_offset;
  const unsigned int bot_iterate_end = bot_end - bot_offset;
  const unsigned int iterate_end = (top_iterate_end > bot_iterate_end) ? \
    bot_iterate_end : top_iterate_end;
  
  for (i = 1; i < iterate_end; i++) {

#if 0
    printf("matches: %d, first_was_match: %d, i: %d, top: %d, bot: %d\n",
	   matches, first_was_match, i,
	   i + top_offset - 1, i + bot_offset - 1);

    std::string s1(top_seq);
    std::string s2(bot_seq);

    s1 = s1.substr(i + top_offset- 1, windowsize);
    s2 = s2.substr(i + bot_offset - 1, windowsize);

    unsigned int test_matches = 0;

    printf("%s\n", s1.c_str());
    printf("%s\n", s2.c_str());
    for (unsigned int k = 0; k < windowsize; k++) {
      if (MATCH(s1, k, s2, k)) {
	test_matches++;
	printf("*");
      } else {
	printf(" ");
      }
    }
    printf("\n");

    assert (test_matches == matches);

#endif

    if (first_was_match) {
      assert(matches >= 1);
      matches--;
    }

    if (MATCH(top_seq, i + top_offset + window_end,
	      bot_seq, i + bot_offset + window_end)) {
      matches++;
    }

    if (matches >= min_matches) {
      match = new Match(top_offset + i, bot_offset + i, windowsize, matches,1);
      if (orientation == -1) { match->reverse_bottom(bot_len); }

      cmp.add_match(match);
    }

    first_was_match = MATCH(top_seq, i + top_offset, bot_seq, i + bot_offset);
  }

  return;
}

void paircomp::_check_parameters(std::string top_seq, std::string bot_seq,
				 unsigned int windowsize, float threshold)
{
  if (top_seq.length() < windowsize ||
      bot_seq.length() < windowsize) {
    throw paircomp_exception("windowsize must be greater than sequence length(s)");
  }

  if (windowsize < 0) {
    throw paircomp_exception("windowsize must be >= 0");
  }
  if (windowsize > 255) {
    throw paircomp_exception("windowsize must be <= 255");
  }

  if (threshold < 0.0 || threshold > 1.0) {
    throw paircomp_exception("threshold must be >= 0 and <= 1.");
  } 
}
