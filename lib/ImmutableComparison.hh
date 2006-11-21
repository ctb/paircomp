#ifndef IMMUTABLE_COMPARISON_HH
#define IMMUTABLE_COMPARISON_HH

//
// Definition of ImmutableComparison class, allowing the manipulation
// of Comparison objects.
//
// None of the ImmutableComparison methods modify the object itself,
// but rather return a new ImmutableComparison.  Use a
// MutableComparison if you need to construct or load a Comparison.
// (You *must* use a MutableComparison to construct a new
// ImmutableComparison.)
//
// This is done so that per-thread object locking need not be done
// for the Python extension.  Well, and it's convenient.
//

#include <string>

#include "Comparison.hh"
#include "MutableComparison.hh"

namespace paircomp {

  // a structure to contain a (fixed) set of matches.
  typedef struct {
    unsigned int num;
    Match * block;
  } _MatchContainer;

  typedef std::map<unsigned int, _MatchContainer* >::const_iterator \
                matches_iterator;

  class ImmutableComparison : public Comparison {
  protected:
    std::map<unsigned int, _MatchContainer *> _matches;

    // utility fn
    void _filter_forward(const ImmutableComparison &bc, bool rev_bc,
			 const ImmutableComparison &ac, bool rev_ac,
			 MutableComparison * new_ab,
			 MutableComparison * new_bc,
			 MutableComparison * new_ac) const;

    // utility fn
    void _build_forward(const ImmutableComparison &bc,
			const char * top, const char * bot,
			const char * bot_rev,
			unsigned int ithreshold,
			MutableComparison &new_ac) const;

    // Protected constructor: create an empty ImmutableComparison.
    ImmutableComparison(unsigned int top, unsigned int bot, unsigned int w)
      : Comparison(top, bot, w) { }
  public:
    // Constructor: takes a MutableComparison only.
    ImmutableComparison(MutableComparison * cmp);

    ~ImmutableComparison();

    // get the _MatchContainer at the given pos, if it exists.
    const _MatchContainer * get_matches(unsigned int pos) const {
      matches_iterator iter = _matches.find(pos);

      iter = _matches.find(pos);
      if (iter != _matches.end()) {
	return iter->second;
      }
      return NULL;
    }

    // Filter the comparison at a given threshold.
    ImmutableComparison * filter_by_threshold(float threshold) const;

    // Filter the comparison for a particular orientation.
    ImmutableComparison * filter_by_orientation(bool forward, bool reverse)
      const;

    // Reverse top or bottom coordinates.
    ImmutableComparison * reverse_top_matches() const;
    ImmutableComparison * reverse_bot_matches() const;

    // Does this comparison contain this match?
    bool contains_match(const Match &m) const;
    
    // Convert the seqcomp into a set of matches between
    // individual bases.
    ImmutableComparison * isolate_matching_bases(std::string top_seq,
						 std::string bot_seq) const;

    // Invert top/bottom.
    ImmutableComparison * invert() const;

    // Subtract out all of the given matches.  Does nothing for
    // Matches that are in 'other' but not 'this', i.e. only the
    // intersection is subtracted.
    ImmutableComparison * subtract(const ImmutableComparison &other) const;

    // Some simple set operations.
    bool is_empty() const;
    bool contains(const ImmutableComparison &other) const;
    bool equals(const ImmutableComparison &other) const;

    ImmutableComparison * intersect(const ImmutableComparison &other) const;
    
    // Output stuff.
    void save_as_seqcomp(char * filename) const;
    void save_as_paircomp(char * filename) const;

    // Three-way stuff.
    ImmutableComparison * build_transitive(const ImmutableComparison &bc,
					   const std::string a_seq,
					   const std::string c_seq,
					   float threshold) const;

    void filter_transitively(const ImmutableComparison &bc,
			     const ImmutableComparison &ac,
			     ImmutableComparison **new_ab,
			     ImmutableComparison **new_bc,
			     ImmutableComparison **new_ac) const;

  };

  class ComparisonFileException : public paircomp_exception {
  public:
    ComparisonFileException(const std::string & m)
      : paircomp_exception(m) {};
  };

};

#endif // IMMUTABLE_COMPARISON_HH
