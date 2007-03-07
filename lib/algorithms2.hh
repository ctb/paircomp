#ifndef ALGORITHMS2_HH
#define ALGORITHMS2_HH

/*!
*   \file algorithms2.hh
*/
#include "algorithms.hh"

/*!
*   \namespace paircomp
*/
namespace paircomp {

  /*! \fn ImmutableComparison * hashed_n_comparison(const std::string seq1,const std::string seq2,unsigned int windowsize,float threshold)
  *   \brief Cuts the top sequence into pieces and finds these pieces in the bottom sequence.
  *   \param seq1 A string object containing the "top" sequence.
  *   \param seq2 A string object containing the "bottom" sequence.
  *   \param windowsize An integer for the scanning window size.
  *   \param threshold A percentage threshold.
  *   \return An ImmutableComparison object with every matching piece.
  */
  ImmutableComparison * hashed_n_comparison(const std::string seq1,
					    const std::string seq2,
					    unsigned int windowsize,
					    float threshold);
}

#endif // ALGORITHMS2_HH
