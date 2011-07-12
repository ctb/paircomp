#ifndef IMMUTABLE_COMPARISON_HH
#define IMMUTABLE_COMPARISON_HH
/*!
*   \file ImmutableComparison.hh
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "Comparison.hh"
#include "MutableComparison.hh"
/*!
*   \namespace paircomp
*/
namespace paircomp {

 
  /*!
  *   \struct _MatchContainer
  *   \brief a structure to contain a (fixed) set of matches.
  */
  typedef struct{
    unsigned int num; /*!< Number of matches in this struct.*/ 
    Match * block; /*!< A pointer to an array of Match objects.*/
  } _MatchContainer;
  
  
  /*!
  *   \var typedef std::map<unsigned int, _MatchContainer* >::const_iterator matches_iterator
  *   \brief This is a mapping of _MatchContainer data structs 
  */
  typedef std::map<unsigned int, _MatchContainer* >::const_iterator matches_iterator;



 /*!
 *   \class ImmutableComparison
 *   \brief Definition of ImmutableComparison class, allowing the manipulation of Comparison objects.
 * 
 *   
 *   None of the ImmutableComparison methods modify the object itself,
 *   but rather return a new ImmutableComparison.  Use a
 *   MutableComparison if you need to construct or load a Comparison.
 *   (You *must* use a MutableComparison to construct a new
 *   ImmutableComparison.)
 *
 *  This is done so that per-thread object locking need not be done
 *  for the Python extension.  Well, and it's convenient.
 *  \sa Comparison
 */
class ImmutableComparison : public Comparison {

  protected:
    /*!
	*   \var std::map<unsigned int, _MatchContainer *> _matches
	*   \brief A mapping of _MatchContainer data structs which is indexable like an array.
	*/
    std::map<unsigned int, _MatchContainer *> _matches;


	/*!
	*   \fn void _filter_forward(const ImmutableComparison &bc, bool rev_bc,
			 const ImmutableComparison &ac, bool rev_ac,
			 MutableComparison * new_ab,
			 MutableComparison * new_bc,
			 MutableComparison * new_ac) const
	*   \brief Function "filters" by looking for transitive "forward" matches between ImmutableComparison objects.
	*   \param bc The first ImmutableComparison object with comparison data for an analysis group.
	*   \param rev_bc A boolean variable to inidcate whether to check the reverse compliment of the bottom seqience of bc.
	*   \param ac The second ImmutableComparison object with comparison data for an analysis group.
	*   \param rev_ac A boolean variable to indicate whether to check the reverse compliment of the bottom seqience of ac.
	*   \param new_ab A pointer to a MutableComparison object which stores the transitive matches for the ab analysis group.
	*   \param new_bc A pointer to a MutableComparison object which stores the transitive matches for the bc analysis group.
	*   \param new_ac A pointer to a MutableComparison object which stores the transitive matches for the ac analysis group.
	*   
	*/
    void _filter_forward(const ImmutableComparison &bc, bool rev_bc,
			 const ImmutableComparison &ac, bool rev_ac,
			 MutableComparison * new_ab,
			 MutableComparison * new_bc,
			 MutableComparison * new_ac) const;


	
	/*!
	*   \fn void _build_forward(const ImmutableComparison &bc,
			const char * top, const char * bot,
			const char * bot_rev,
			unsigned int ithreshold,
			MutableComparison &new_ac) const
	*   \brief Builds a list of forward matches and places them into the MutableComparison parameter new_ac.
	*   \param bc A reference to an ImmutableComparison object, that represents Matches for an analysis group.
	*   \param top A character array representing the "top" sequence.
	*   \param bot A character array representing the "bottom" sequence.
	*   \param bot_rev A boolean indicating whether to reverse the "bottom" sequence.
	*   \param ithreshold An integer threshold 
	*   \param new_ac A reference to an ImmutableComparison object that represents the Matches for an analysis group.
	*
	*  Function builds a list of ImmutableComparison objects using forward matches from a->b
	*  and places the Match objects which have Match::_n above the value of ithreshold,
	*  into the MutableComparison object new_ac.
	*/
    void _build_forward(const ImmutableComparison &bc,
			const char * top, const char * bot,
			const char * bot_rev,
			unsigned int ithreshold,
			MutableComparison &new_ac) const;


	/*!
	*   \param top An integer indicating the desired length of the top sequence.
	*   \param bot An integer indicating the desired length of the bottom sequence.
	*   \param w An integer indicating the windowsize for this object.
	*  
	*   Protected constructor: create an empty ImmutableComparison.
	*/
    ImmutableComparison(unsigned int top, unsigned int bot, unsigned int w)
      : Comparison(top, bot, w) { }
  public:
	
	/*!
	*   \param cmp A pointer to a MutableComparison to be transformed into an ImmutableComparison.
	*   
	*   Constructor: takes a MutableComparison & transforms it.
	*/
    ImmutableComparison(MutableComparison * cmp);

    ~ImmutableComparison();


	/*!
	*   \fn const _MatchContainer * get_matches(unsigned int pos) const
	*   \brief Get the _MatchContainer at the given pos, if it exists.
	*   \param pos An integer value indicating the position of the _MatchContainer object in _matches to return.
	*   \return The _MatchContainer object at position pos, false if there is nothing at that position.
	*
	*/
    const _MatchContainer * get_matches(unsigned int pos) const {
      matches_iterator iter = _matches.find(pos);

      iter = _matches.find(pos);
      if (iter != _matches.end()) {
	     return iter->second;
      }
      return NULL;
    }


	/*!
	*   \fn ImmutableComparison * filter_by_threshold(float threshold) const
	*   \brief Filter the comparison at a given threshold.
	*   \param threshold A percentage threshold a match should pass (is rounded up).
	*   \return An ImmutableComparison object with all of the matches equal/above the threshold.
	*
	*   Function goes through all matches in _matches. Matches above the threshold 
	*   are placed into a new ImmutableComparison object.
	*/
    ImmutableComparison * filter_by_threshold(float threshold) const;

	/*!
	*   \fn ImmutableComparison * filter_by_orientation(bool fwd, bool rev) const
	*   \brief Filter the comparison to contain only matches in the given orientation.
	*   \param fwd Only retain forward matches
	*   \param rev Only retain reverse matches
	*   \return An ImmutableComparison object with all of the matches in the right orientation.
	*
	*   Function goes through all matches in _matches. Matches with the correct orientation(s)
	*   are placed into a new ImmutableComparison object.
	*/
    ImmutableComparison * filter_by_orientation(bool fwd, bool rev) const;

	/*!
	*   \fn ImmutableComparison * reverse_top_matches() const
	*   \brief Reverse top coordinates.
	*   \return Returns a copy of ImmutableComparison with the top "reversed."
	*
	*   For each Match in _matches, the top sequence is reversed.
	*/
    ImmutableComparison * reverse_top_matches() const;
	
	/*!
	*   \fn ImmutableComparison * reverse_bot_matches() const
	*   \brief Reverse bottom coordinates.
	*   \return Returns a copy of ImmutableComparison with the bottom "reversed."
	*
	*   For each Match in _matches, the bottom sequence is reversed.
	*/
    ImmutableComparison * reverse_bot_matches() const;


	/*!
	*   \fn bool contains_match(const Match &m) const
	*   \brief Does this comparison contain this match?
	*   \param m A reference to a Match object to check for.
	*   \return True if m is contained in _matches, false otherwise.
	*/
    bool contains_match(const Match &m) const;
    
  
	/*!
	*   \fn ImmutableComparison * isolate_matching_bases(std::string top_seq,
						 std::string bot_seq) const
	*   \brief Convert the seqcomp into a set of matches between individual bases.
	*   \param top_seq A string object holding the top sequence.
	*   \param bot_seq A string object holding the bottom sequence.
	*   \return An ImmutableComparison object with has a Match for each matching bp.
	*   \exception paircomp_exception lengths do not match
	*
	*   Convert the comparison into a set of 1-bp matches between actually-matching
    *   bases, given the sequences to compare. 
	*/
    ImmutableComparison * isolate_matching_bases(std::string top_seq,
						 std::string bot_seq) const;


	/*!
	*   \fn ImmutableComparison * invert() const
	*   \brief Inverts the Match objects in _matches and returns them in a new ImmutableComparison. 
	*   \return A copy of this ImmutableComparison with Top/Bot inverted.
	*   \sa Match::invert
	*/
    ImmutableComparison * invert() const;

  
	/*!
	*   \fn ImmutableComparison * subtract(const ImmutableComparison &other) const
	*   \brief Subtract out all of the given matches. 
	*   \param other A reference to another ImmutableComparison object.
	*   \return A pointer to an ImmutableComparison with matches that are in this object but not other. 
	*
	*   Does nothing for Matches that are in 'other' but not 'this', i.e. only the 
	*   intersection is subtracted.
	*/
    ImmutableComparison * subtract(const ImmutableComparison &other) const;


	/*!
	*   \fn bool is_empty() const
	*   \brief Checks if there are any matches in this object.
	*   \return True if there is at least one match, false otherwise.
	*/
    bool is_empty() const;
	
	/*!
	*   \fn bool contains(const ImmutableComparison &other) const
	*   \brief Checks to see if parameter other is contained in this ImmutableComparison.
	*   \param other A reference to another ImmutableComparison object.
	*   \return True if the matches in "other" are present within this object, false otherwise.
	*   \sa subtract
	*   
	*   Function Subtracts the matches in this object from other and checks to 
	*   see if the resulting ImmutableComparison is empty.
	*/
    bool contains(const ImmutableComparison &other) const;
	
	
	/*!
	*   \fn bool equals(const ImmutableComparison &other) const
	*   \brief Checks to see if parameter other is equal to this ImmutableComparison.
	*   \param other A reference to another ImmutableComparison.
	*   \return True if "other" is equal to this object, false otherwise.
	*   \sa contains
	*
	*   Calls function contains to check for equality. 
	*/
    bool equals(const ImmutableComparison &other) const;



    /*!
	*   \fn ImmutableComparison * intersect(const ImmutableComparison &other) const
	*   \brief Returns all intersecting matches.
	*   \param other A reference to another ImmutableComparison.
	*   \return An ImmutableComparison containing all interecting Match objects.
	*
	*   Function saves all Matches that are in both this ImmutableComparison
	*   and "other."
	*/
    ImmutableComparison * intersect(const ImmutableComparison &other) const;
    

	/*!
	*   \fn void save_as_seqcomp(char * filename) const
	*   \brief Output a seqcomp formated file.
	*   \param filename A character array with a filename to write the output to.
	*
	*   Save to a given file in Tristan-style b3.5 output.
	*/
    void save_as_seqcomp(char * filename) const;
	
	/*!
	*   \fn void save_as_paircomp(char * filename) const
	*   \brief Output a paircomp formated file.
	*   \param filename A character array with a filename to write the output to.
	*
	*   This outputs the raw data for each match row by row into a text file.
	*      >top_pos	bot_pos	n_matching	orientation
	*/
    void save_as_paircomp(char * filename) const;



    // Three-way stuff.
	/*!
	*   \fn ImmutableComparison * build_transitive(const ImmutableComparison &bc,
					   const std::string a_seq,
					   const std::string c_seq,
					   float threshold) const
	*   \brief Builds a transitive ImmutableComparison from two sequences and another ImmutableComparison. 
	*   \param bc An ImmutableComparison that contains the matches for the "bc" Analysis group to use to build the transitive ImmutableComparison.
	*   \param a_seq A string objects containing the "a" sequence to building the transitive ImmutableComparison.
	*   \param c_seq A string objects containing the "b" sequence to building the transitive ImmutableComparison.
	*   \param threshold The threshold for the minimum number of matches per _windowsize to use when building.
	*   \return A pointer to an ImmutableComparison with The transitive Match objects between a_seq, b_seq, and bc.
	*   \sa _build_forward
	*
	*   Function uses _build_forward to generate matches between a_seq, b_seq, bc and the reversal of bc.
	*/
    ImmutableComparison * build_transitive(const ImmutableComparison &bc,
					   const std::string a_seq,
					   const std::string c_seq,
					   float threshold) const;
    
	/*!
	*   \fn void filter_transitively(const ImmutableComparison &bc,
			     const ImmutableComparison &ac,
			     ImmutableComparison **new_ab,
			     ImmutableComparison **new_bc,
			     ImmutableComparison **new_ac) const
    *   \brief Function filters transitive matches from bc and ac and creates new ImmutableComparison objects new_ab, new_bc, and new_ac.
	*   \param bc An ImmutableComparison with matches for a bc AnalysisGroup.
	*   \param ac An ImmutableComparison with matches for an ac AnalysisGroup.
	*   \param new_ab A pointer to a pointer of a new ImmutableComparison for storing ab Matches filtered from bc and ac.
	*   \param new_bc A pointer to a pointer of a new ImmutableComparison for storing bc Matches filtered from bc and ac.
	*   \param new_ac A pointer to a pointer of a new ImmutableComparison for storing ac Matches filtered from bc and ac.
	*   \sa _filter_forward
	*
	*   Function uses _filter_forward to filter through bc and ac to obtain transitive matches
	*   and places them in new_ab, new_bc, and new_ac respectively. The args new_ab, new_bc, and new_ac
	*   are pointers to pointers, so the address of the pointer to the new ImmutableComparison
	*   must be passed to this function.
	*
	*/
    void filter_transitively(const ImmutableComparison &bc,
			     const ImmutableComparison &ac,
			     ImmutableComparison **new_ab,
			     ImmutableComparison **new_bc,
			     ImmutableComparison **new_ac) const;

  };
  
  
  /*!
  *   \class ComparisonFileException
  *   \brief A class for throwing exceptions with Comparisons.
  *   \sa paircomp_exception
  */
  class ComparisonFileException : public paircomp_exception {
  public:
  
    /*!
	*   \param m A string object for an error message.
	*/
    ComparisonFileException(const std::string & m)
      : paircomp_exception(m) {};
  };

};

#endif // IMMUTABLE_COMPARISON_HH
