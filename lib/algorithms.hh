/*!
*   \file algorithms.hh
*/
#ifndef ALGORITHMS_HH
#define ALGORITHMS_HH

#include <string>
#include "MutableComparison.hh"
#include "ImmutableComparison.hh"
/*!
*   \namespace paircomp
*/
namespace paircomp {

  /*! \fn ImmutableComparison * simple_nxn_comparison(const std::string seq1,const std::string seq2,unsigned int windowsize,float threshold)
  *   \brief performs a simple NxMxW comparison, function performs an N squared forward comparison, then an N squared reverse comparison.
  *   \param seq1 The top DNA/RNA/Protein sequence to be used for a comparison.
  *   \param seq2 The bottom DNA/RNA/Protein sequence to be used for comparison.
  *   \param windowsize The number of characters to compare in the scanning window.
  *   \param threshold Number of matches that should be present in a particular scanning window.
  *   \return An ImmutableComparison object containing the match data from the simple comparison.
  *   \sa ImmutableComparison
  */
  ImmutableComparison * simple_nxn_comparison(const std::string seq1,
					      const std::string seq2,
					      unsigned int windowsize,
					      float threshold);
						  
  
  /*! \fn ImmutableComparison * rolling_nxn_comparison(const std::string seq1,const std::string seq2,unsigned int windowsize,float threshold)
  *   \brief Function performs a linear comparison across two sequences.
  *   \param seq1 The top DNA/RNA/Protein sequence to be used for a comparison.
  *   \param seq2 The bottom DNA/RNA/Protein sequence to be used for comparison.
  *   \param windowsize The number of characters to compare in the scanning window.
  *   \param threshold Number of matches that should be present in a particular scanning window.
  *   \return An ImmutableComparison object containing the match data from the rolling comparison
  *   \sa ImmutableComparison
  */
  ImmutableComparison * rolling_nxn_comparison(const std::string seq1,
					       const std::string seq2,
					       unsigned int windowsize,
					       float threshold);

  //
  // Some useful utility functions.
  //
  
  /*! \fn void void _preprocess(std::string &seq)
  *   \brief Converts a sequence into all uppercase.
  *   \param seq starting address to a DNA/RNA/Protein sequence to be converted to uppercase.
  */
  // transform sequence to uppercase.
  void _preprocess(std::string &seq);
  
  
  
  
  /*! \fn std::string _reverse_complement(const std::string seq)
  *   \brief Obtains the reverse compliment of a given input sequence.
  *   \param seq The address of a string with a DNA/RNA/Protein sequence.
  *   \return A reverse compliment copy of the original string. 
  */
  std::string _reverse_complement(const std::string seq);




  /*! \fn bool _match_window(const char * top, unsigned int start_top,const char * bot, unsigned int start_bot,const unsigned int windowsize,const unsigned int threshold,unsigned int &matches)
  *   \brief Funtion scans a "window" between two sequences and counts for a number of matches present between the two.
  *   \param top A pointer to a character array representing a top DNA/RNA/Protein sequence for comparison.
  *   \param start_top An integer indicating a start index to start from in the top sequence.
  *   \param bot A pointer to a character array representing a bottom DNA/RNA/Protein sequence for comparison
  *   \param start_bot An integer indicating a start index to start from in the bottom sequence.
  *   \param windowsize An integer indicating a scanning window size for comparing the sequences.
  *   \param threshold An integer which indicates the number of matches that should be present for a window to match another.
  *   \param matches The number of matches present in a window comparison.
  *   \return Returns true if the number of matches is greater than or equal to the threshold. False otherwise.
  */
  bool _match_window(const char * top, unsigned int start_top,
		     const char * bot, unsigned int start_bot,
		     const unsigned int windowsize,
		     const unsigned int threshold,
		     unsigned int &matches);
  
  
  
  
  /*! \fn void _compare_stretch(MutableComparison &cmp,unsigned int windowsize,unsigned int min_matches,const char * top_seq, unsigned int top_end,unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const int orientation)
  *   \brief Function compares a window between two sequences from an of set to an end point and find matches, matches over the min_matches are loaded into the MutableComparison which is passed by reference.
  *   \param cmp A pointer to a MutableComparison object for adding match data to.
  *   \param windowsize An integer indicating the scanning window size.
  *   \param min_matches The minumum number of matches which must be attained before adding a match object to the MutableComparison object (cmp)
  *   \param top_seq A character array representing the top DNA/RNA/Protein sequence.
  *   \param top_end An integer index into the top sequence indicating a stop point for scanning.
  *   \param top_offset An integer offset to the start of the top sequence.
  *   \param bot_seq A character array representing the bottom DNA/RNA/Protein sequence.
  *   \param bot_end An integer index into the bottom sequence indicating a stop point for scanning.
  *   \param bot_offset An integer offset to the start of the bottom sequence.
  *   \param orientation A flag indicating whether to reverse the bottom sequence of a match.  -1 to flag for reverse, no reverse otherwise.
  *   \sa Match
  *
  *   Compares two stretches between two sequences and records the result of a window scan if the number
  *   of matches are above min_matches. If orientation == -1 then the bottom sequence of the match to be returned
  *   via reference is reversed. 
  */
  void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const int orientation);



  /*! \fn void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const char * bot_seq_rev)
  *   \brief Function compares a window between two sequences, it returns a referenced MutableComparison object which contains matches with a forward and reversed bottom sequence.
  *   \param cmp A pointer to a MutableComparison object for adding match data to.
  *   \param windowsize An integer indicating the scanning window size.
  *   \param min_matches The minumum number of matches which must be attained before adding a match object to the MutableComparison object (cmp)
  *   \param top_seq A character array representing the top DNA/RNA/Protein sequence.
  *   \param top_end An integer index into the top sequence indicating a stop point for scanning.
  *   \param top_offset An integer offset to the start of the top sequence.
  *   \param bot_seq A character array representing the bottom DNA/RNA/Protein sequence.
  *   \param bot_end An integer index into the bottom sequence indicating a stop point for scanning.
  *   \param bot_offset An integer offset to the start of the bottom sequence.
  *   \param bot_seq_rev A character array representing a REVERSED bottom DNA/RNA/Protein sequence.
  *   \sa void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const int orientation)
  *   
  *   Calls the like named function and performs a stretch comparison for both orientations.
  *   
  */
  void _compare_stretch(MutableComparison &cmp,
			unsigned int windowsize,
			unsigned int min_matches,
			const char * top_seq, unsigned int top_end,
			unsigned int top_offset,
			const char * bot_seq, unsigned int bot_end,
			unsigned int bot_offset,
			const char * bot_seq_rev);
  
  
  /*! \fn unsigned int _n_matching(const char * a, unsigned int a_start,
			   const char * b, unsigned int b_start,
			   unsigned int w)
  *   \brief Performs A simple linear comparison, counts the number of matches between the windows of a and b.
  *   \param a A character array representing the first DNA/RNA/Protein sequence.
  *   \param a_start The starting index of character array a for linear comparison.
  *   \param b A character array presenting the second DNA/RNA/Protein sequence.
  *   \param b_start The starting index of character array b for linear comparison.
  *   \param w The scanning window size.
  *   \return The number of character matches between sequence a and b.
  *
  *   Function compares a window of 'w' number of characters between two character 
  *   arrays and returns the number of matches.
  */
  unsigned int _n_matching(const char * a, unsigned int a_start,
			   const char * b, unsigned int b_start,
			   unsigned int w);
  
  
  /*! \fn void _check_parameters(std::string top_seq, std::string bot_seq,
			 unsigned int windowsize, float threshold)
  *   \brief Function checks parameters given to a rolling and nxn comparison.
  *   \param top_seq The top DNA/RNA/Protein sequence.
  *   \param bot_seq The bottom DNA/RNA/Protein sequence.
  *   \param windowsize An integer indicating the scanning window size.
  *   \param threshold A float value indicating the threshold percentage.
  *   \exception paircomp_exception windowsize must be greater than sequence length(s)
  *   \exception paircomp_exception windowsize must be >= 0
  *   \exception paircomp_exception windowsize must be <= 255
  *   \exception paircomp_exception threshold must be >= 0 and <= 1.
  *
  *   Function throws an exception if parameters are invalid.
  *
  *
  */
  void _check_parameters(std::string top_seq, std::string bot_seq,
			 unsigned int windowsize, float threshold);
}

#endif // ALGORITHMS_HH
