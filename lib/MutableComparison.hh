#ifndef MUTABLE_COMPARISON_HH
#define MUTABLE_COMPARISON_HH
/*!
*   \file MutableComparison.hh
*/
#include <vector>
#include <map>

#include "Match.hh"
#include "Comparison.hh"

/*!
*   \namespace paircomp
*/
namespace paircomp {
  class ImmutableComparison;

  
  
  /*!
   *   \class MutableComparison
   *
   *   \brief Definition of MutableComparison class, for constructing fixed-width window comparisons between two DNA sequences.
   *   \sa ImmutableComparison
   *   \sa Comparison
   *
   *   MutableComparison is a class especially designed for adding new
   *   matches; ImmutableComparison should be used for manipulation etc.
   *
   *   Use MutableComparison when creating new analyses or loading analyses
   *   from a file.
   */
  class MutableComparison : public Comparison {
    friend class ImmutableComparison;
  protected:
  
    /*!
	*   \var _matches
	*   \brief A mapping of Match object vectors which is indexable like an array.
	*
	*   This data member stores all of the Match data objects for this class.
	*/
    std::map<unsigned int, std::vector<Match*> > _matches;
  public:
  
    /*!
	*   \brief Constructor: calls the Comparison constructor for it's underlying data structure.
	*   \param top The length of the top sequence for this comparison data.
	*   \param bot The length of the bottom sequence for this comparison data.
	*   \param w The window size used for this comparison data.
	*/
    MutableComparison(unsigned int top, unsigned int bot, unsigned int w);
	
	
	/*!
	*
	*/
    ~MutableComparison();


	/*!
	*   \fn void add_match(Match * &match, bool require_unique=true)
	*   \brief Add the match to the given position.  (Consumes the match.)
	*   \param match A reference to a Match to add to the MutableComparison::_matches class member.
	*   \param require_unique Forces the function to only add unique matches when true. Default is set to true.
	*   \exception paircomp_exception match is not unique, as promised!
	*/
    void add_match(Match * &match, bool require_unique=true);

    
	/*!
	*   \fn void parse_paircomp_format(char * buf)
	*   \brief Load from various formats.
	*   \param buf An array of characters read from a file containing text in paircomp format to be parsed.
	*   \exception ComparisonParserException top length as given is incorrect
	*   \exception ComparisonParserException bottom match out of bounds
	*   \exception ComparisonParserException n matches larger than windowsize!
	*/
    void parse_paircomp_format(char * buf);
	
	/*!
	*   \fn void parse_seqcomp_format(const char * buf)
	*   \brief Parses a seq_comp formatted (b3.5 format) array of characters.
	*   \param buf An array of characters read from a file containing text in b3.5 format to be parsed.
	*   \exception ComparisonParserException magic number 'b3.5' not present in file
	*   \exception ComparisonParserException windowsize (in buf) does not match expected windowsize
	*   \exception ComparisonParserException seqcomp lines in buf do not match expected top length.
	*   \exception ComparisonParserException unable to understand orientation
	*   \exception ComparisonParserException bottom match out of bounds
	*   \exception ComparisonParserException n matches larger than windowsize!
    */
    void parse_seqcomp_format(const char * buf);
    
}; //ends namespace paircomp	
  
  /*!
  *   \fn MutableComparison * load_paircomp_file(unsigned int top_len, unsigned int bot_len, 
                    unsigned int windowsize, char * filename)
  *   \brief Function loads a paircomp file into a MutableComparison Data object.
  *   \param top_len The length of the top sequence for this MutableComparison to be loaded.
  *   \param bot_len The length of the top sequence for this MutableComparison to be loaded.
  *   \param windowsize The windowsize for this MutableComparison to be loaded.
  *   \param filename A character array with a filename to open.
  *   \return A pointer to a new MutableComparison with the data parsed from the file.
  *
  */
  MutableComparison * load_paircomp_file(unsigned int top_len, unsigned int bot_len, 
		unsigned int windowsize, char * filename);
  
  
  /*!
  *   \fn MutableComparison * load_seqcomp_file(unsigned int top_len, unsigned int bot_len, 
                   unsigned int windowsize, char * filename)
  *   \brief Function loads a paircomp file into a MutableComparison Data object.
  *   \param top_len The length of the top sequence for this MutableComparison to be loaded.
  *   \param bot_len The length of the top sequence for this MutableComparison to be loaded.
  *   \param windowsize The windowsize for this MutableComparison to be loaded.
  *   \param filename A character array with a filename to open.
  *   \return A pointer to a new MutableComparison with the data parsed from the file.
  *
  */
  MutableComparison * load_seqcomp_file(unsigned int top_len, unsigned int bot_len, 
           unsigned int windowsize, char * filename);

  
  /*!
  *   \class ComparisonParserException
  *   \brief A class for throwing exceptions when Comparison parsing.
  *   \sa paircomp_exception
  */
  class ComparisonParserException : public paircomp_exception {
  public:
  
    /*!
	*   \param m A reference to a string with an error message. 
	*/
    ComparisonParserException(const std::string & m)
      : paircomp_exception(m) { ; };
 //};//fake end to the paircomp namespace for doxygen to parse these functions, just comment this out before compiling, and uncomment the previous namespace end. 

  };
};

#endif // MUTABLE_COMPARISON_HH
