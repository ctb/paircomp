#ifndef MUTABLE_COMPARISON_HH
#define MUTABLE_COMPARISON_HH

//
// Definition of MutableComparison class, for constructing fixed-width
// window comparisons between two DNA sequences.
//
// MutableComparison is a class especially designed for adding new
// matches; ImmutableComparison should be used for manipulation etc.
//
// Use MutableComparison when creating new analyses or loading analyses
// from a file.
//

#include <vector>
#include <map>

#include "Match.hh"
#include "Comparison.hh"

namespace paircomp {
  class ImmutableComparison;

  class MutableComparison : public Comparison {
    friend class ImmutableComparison;
  protected:
    std::map<unsigned int, std::vector<Match*> > _matches;
  public:
    MutableComparison(unsigned int top, unsigned int bot, unsigned int w);
    ~MutableComparison();

    // Add the match to the given position.  (Consumes the match.)
    void add_match(Match * &match, bool require_unique=true);

    // Load from various formats.
    void parse_paircomp_format(char * buf);
    void parse_seqcomp_format(const char * buf);
  };

  MutableComparison * load_paircomp_file(unsigned int top_len,
					 unsigned int bot_len,
					 unsigned int windowsize,
					 char * filename);

  MutableComparison * load_seqcomp_file(unsigned int top_len,
					unsigned int bot_len,
					unsigned int windowsize,
					char * filename);

  class ComparisonParserException : public paircomp_exception {
  public:
    ComparisonParserException(const std::string & m)
      : paircomp_exception(m) { ; };
  };
};

#endif // MUTABLE_COMPARISON_HH
