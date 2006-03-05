#ifndef COMPARISON_HH
#define COMPARISON_HH

//
// Definition of basic Comparison class, containing matches from
// a fixed-width window comparison between two DNA sequences.
//
// This is the base class for MutableComparison and
// ImmutableComparison, which are classes that allow the construction
// & manipulation of comparisons, respectively.
//
// The reason for this division between Immutable and Mutable
// comparison classes is that the underlying representation of the
// classes is optimized for their respective roles.  This may strike
// you as a half-assed reason but it's all I got ;).
//
// (In particular, memory fragmentation is definitely an issue.)
//
// Actually, another great reason is that since ImmutableComparison
// methods never modify the object, no per-thread object locking
// needs to be done for the Python interface.
//

#include <exception>
#include <string>

namespace paircomp {
  class paircomp_exception : public std::exception {
  public:
    const std::string error_msg;
    paircomp_exception(const std::string s) : error_msg(s) {};
    ~paircomp_exception() throw () {};
  };

  class Comparison {
  protected:
    const unsigned int _top_length, _bottom_length, _windowsize;

    // protected constructor: allow only subclasses to initialize
    Comparison(unsigned int top_len,
	       unsigned int bot_len,
	       unsigned int winsize) : _top_length(top_len),
				       _bottom_length(bot_len),
				       _windowsize(winsize) { }
  public:
    unsigned int get_top_length() const { return _top_length; }
    unsigned int get_bottom_length() const { return _bottom_length; }
    unsigned int get_windowsize() const { return _windowsize; }

    virtual ~Comparison() { };
  };
};

#endif // COMPARISON_HH
