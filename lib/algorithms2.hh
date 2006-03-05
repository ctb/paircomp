#ifndef ALGORITHMS2_HH
#define ALGORITHMS2_HH

#include "algorithms.hh"

namespace paircomp {
  ImmutableComparison * hashed_n_comparison(const std::string seq1,
					    const std::string seq2,
					    unsigned int windowsize,
					    float threshold);
}

#endif // ALGORITHMS2_HH
