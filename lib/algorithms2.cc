#include "algorithms2.hh"

/*********************************************************************
 *
 * This code (algorithms2.cc) is modified from the Dot Plot Sequence Comparison
 * Code.  The original author is Shoudan Liang, and the code was
 * subsequently modified by Titus Brown for inclusion in the paircomp
 * subproject of the FamilyJewels project.
 *
 * The Dot Plot Sequence Comparison Code was developed by  the U.S.
 * Government as represented by the Administrator of the National
 * Aeronautics and Space Administration.  No copyright is claimed in the
 * United States on behalf of the U.S. government.
 *
 * The Dot Plot Sequence Comparison Code may be used, copied, and
 * provided to others as part of the FamilyJewels - Comparative sequence
 * analysis for genomic data  software distributed by California Institute
 * of Technology under the GNU General Public License.  This software is
 * provided without any warranty, either express or implied.
 *
 * For more information about  The Dot Plot Sequence Comparison Code
 * contact Shoudan Liang, NASA Ames Research Center, Moffett Field, CA
 * 94035 (E-mail: Shoudan.Liang@nasa.gov, Phone: 650-604-6631)
 *
 *********************************************************************/

using namespace paircomp;

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <iostream>

static void _find_match(const unsigned int top_pos,
			const std::string current_pattern,
			const unsigned int windowsize,
			const unsigned int hard_threshold,
			std::set<std::string> &substring,
			std::map<std::string, std::set<unsigned int> > &word_locations,
			int mismatches_allowed,
			std::string &word,
			MutableComparison &cmp);

ImmutableComparison * paircomp::hashed_n_comparison(const std::string seq1,
						    const std::string seq2,
						    unsigned int windowsize,
						    float threshold)
{
  unsigned int hard_threshold;
  std::set<std::string> substring;

  // each set contains the locations of the substring
  std::map<std::string, std::set<unsigned int> > word_locations;

  // transform the strings appropriately: convert to uppercase, etc.
  std::string top_seq = seq1;
  std::string bot_seq = seq2;

  _preprocess(top_seq);
  _preprocess(bot_seq);
  
  _check_parameters(top_seq, bot_seq, windowsize, threshold);

  hard_threshold = (unsigned int) (threshold * (float)windowsize + 0.5);

  // build an associate array from bot_seq
  for (unsigned int i2=0; i2+windowsize <= bot_seq.size(); i2++) {
    // cut a word of length windowsize from bot_seq starting from position i2
    std::string word = bot_seq.substr(i2,windowsize);
    word_locations[word].insert(i2);

    // create the substring tree
    for (unsigned int i=1; i < windowsize; i++) {
      std::string sub_word = word.substr(0, i);
      substring.insert(sub_word);
    }
  }

  MutableComparison cmp(top_seq.length(),
			bot_seq.length(),
			windowsize);

  // cut top_seq into short words and find all locations of the word in bot_seq
  for (unsigned int i1=0; i1 + windowsize <= top_seq.size(); i1++) {
    std::string word;

    _find_match(i1, top_seq.substr(i1, windowsize),
		windowsize, hard_threshold,
		substring, word_locations,
		windowsize-hard_threshold, word, cmp);

  }

  return new ImmutableComparison(&cmp);
}

static void _find_match(const unsigned int top_pos,
			const std::string current_pattern,
			const unsigned int windowsize,
			const unsigned int hard_threshold,
			std::set<std::string> &substring,
			std::map<std::string, std::set<unsigned int> > &word_locations,
			int mismatches_allowed, std::string &word,
			MutableComparison &cmp) {
  static char nucl[]={'A','C','G','T', 'N'};
  if (word.size() == windowsize) {		// done! output the results
    std::set<unsigned int>::iterator it = word_locations[word].begin();

    for (; it != word_locations[word].end(); it++) {
      Match * m = new Match(top_pos, *it, windowsize, hard_threshold + mismatches_allowed, 1);
      cmp.add_match(m);
    }

    return;
  }

  char c = current_pattern[word.size()];
  word += c;

  // full-length word is not in substring
  if (substring.count(word) || word.size() == windowsize) {
    if ('N'==c) {
      if (mismatches_allowed > 0) {
	// 'N' doesn't match anything
	_find_match(top_pos, current_pattern, windowsize, hard_threshold, substring, word_locations, mismatches_allowed-1, word, cmp);
      }
    } else {
      _find_match(top_pos, current_pattern, windowsize, hard_threshold, substring, word_locations, mismatches_allowed, word, cmp);
    }
  }
  if (mismatches_allowed-- > 0) {
    for (int t=0; t<5; t++) {
      if (nucl[t] == c)		// the case of no mismatch, done already 
	continue;
      word[word.size()-1]=nucl[t];

      // full-length word is not in substring
      if (substring.count(word) || word.size() == windowsize)
	_find_match(top_pos, current_pattern, windowsize, hard_threshold, substring, word_locations, mismatches_allowed, word, cmp);
    }
  }
  // erase the last charactor of word, restoring word to the
  // original input, this also reduce the size of word by 1
  word.erase(word.end()-1);

  return;
}
