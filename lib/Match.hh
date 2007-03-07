/*!
*   \file Match.hh
*/

// Class definition for fixed-width window match.
//
// See README.txt for license and copyright information.

#ifndef MATCH_HH		// guard against multiple inclusion
#define MATCH_HH

#include <assert.h>

/*!
*   \class Match
*   \brief Class holds data and operations relevant to a Match between to sequences.
*
*   Note, this class does not contain the actual sequences as data, only the 
*   start positions, length, number of matches and orientation.
*/
class Match {
protected:
  /*!
  *   \var unsigned int _tpos
  *   \brief Index into the start of a "top" sequence where a match was found.
  */
  unsigned int _tpos; 
  
  /*!
  *   \var unsigned int _bpos
  *   \brief Index into the start of a "bottom" sequence where a match was found.
  */
  unsigned int _bpos;
  
  /*!
  *   \var unsigned char _l
  *   \brief The length of this "match."
  */
  unsigned char _l;
  
  /*!
  *   \var unsigned char _n
  *   \brief The number of character matches.
  */
  unsigned char _n;
  
  /*!
  *   \var char _orientation
  *   \brief The orientation of the bottom sequence.
  *
  *   1 is for forward oriented.
  *   -1 is for reverse oriented.
  */
  char _orientation;
public:

  /*!
  *   \param top_pos Index to the top sequence for the start of a match.
  *   \param bot_pos Index to the bottom sequence for the start of a match.
  *   \param length The length of this match.
  *   \param n_matching Number of matches between two sequences.
  *   \param orientation the orientation (forward or reversed) of the sequences being matched.
  */
  Match(unsigned int top_pos, unsigned int bot_pos,
	unsigned int length, unsigned int n_matching,
	int orientation) {
    assert(length < 256);

    assert(n_matching < 256);
    assert(n_matching <= length);

    assert(orientation == -1 || orientation == 1);

    _tpos = top_pos;
    _bpos = bot_pos;
    _l = length;
    _n = n_matching;
    _orientation = orientation;
  };
  

  /*!
  *   \fn unsigned int get_top_pos() const
  *   \brief Returns the top position index.
  *   \return An unsigned int, for the top position.
  *   \sa _tpos
  */
  unsigned int get_top_pos() const { return _tpos; };
  
  /*!
  *   \fn unsigned int get_bot_pos() const
  *   \brief Returns the bottom position index.
  *   \return An unsigned int, for the bottom position.
  *   \sa _bpos
  */
  unsigned int get_bot_pos() const { return _bpos; };
  
  /*!
  *   \fn unsigned int get_length() const
  *   \brief Returns the length of this match.
  *   \return An unsigned int for the length.
  *   \sa _l
  *  
  *   Note, actual data member is a char, it's returned as an unsigned int.
  */
  unsigned int get_length() const { return (unsigned int) _l; };
  
  /*!
  *   \fn unsigned int get_n_matching() const
  *   \brief Returns the number of matches.
  *   \return An unsigned int for the number of matches. 
  *   \sa _n
  *
  *   Note, actual data member is a char, it's returned as an unsigned int.  
  */
  unsigned int get_n_matching() const { return (unsigned int) _n; };
  
  /*!
  *   \fn int get_orientation() const
  *   \brief Returns the match orientation.
  *   \return An int with a value of either 1 or -1.
  *   \sa _orientation
  */
  int get_orientation() const { return (int) _orientation; };

  /*!
  *   \fn Match * copy() const
  *   \brief Copies all of the data from this Match and returns a copy of it. 
  *   \return A pointer to a copy of this Match object.
  */
  Match * copy() const {
    return new Match(_tpos, _bpos, _l, _n, _orientation);
  }


  /*!
  *   \fn void reverse_top(unsigned int top_len, unsigned int bot_len)
  *   \brief Function "reverses" the coordinates of the "top" sequence of this Match.
  *   \param top_len The length of the "top" sequence of this Match.
  *   \param bot_len The length of the "bottom" sequence of this Match.
  */
  void reverse_top(unsigned int top_len, unsigned int bot_len) {
    if (_orientation == -1) {	// good gawd
      _bpos = _bpos - (_l - 1);
    } else {
      _bpos = _bpos + (_l - 1);
    }

    _tpos = top_len - _tpos - _l;
    _orientation = -1 * _orientation;
  }

  /*!
  *   \fn void reverse_bottom(unsigned int bot_len)
  *   \brief Function "reverses" the coordinates of the "bottom" sequence of this Match. 
  *   \param bot_len The length of the "bottom" sequence of this Match.
  */
  void reverse_bottom(unsigned int bot_len) {
    _bpos = bot_len - (_bpos + 1);
    _orientation = -1 * _orientation;
  }


  /*!
  *   \fn void invert(unsigned int top_len, unsigned int bot_len)
  *   \brief Inverts this Match by swaping values _bpos and _tpos
  *   \param top_len
  *   \param bot_len
  *
  *   Function doesn't seem to use the parameters. :/
  *   Function safely swaps values with respect to orientation.
  */
  void invert(unsigned int top_len, unsigned int bot_len) {
    if (_orientation == -1) {
      // reverse bottom
      unsigned int new_tpos = _bpos - (_l - 1);

      // reverse top
      unsigned int new_bpos = _tpos + (_l - 1);

      // invert
      _tpos = new_tpos;
      _bpos = new_bpos;
    } else {
      unsigned int new_tpos = _bpos;
      unsigned int new_bpos = _tpos;

      _tpos = new_tpos;
      _bpos = new_bpos;
    }
  }


  /*!
  *   \fn bool equals(const Match * other) const
  *   \brief Checks to see if parameter 'other' is equal to this object.
  *   \param other A pointer to another Match object.
  *   \return True if 'other' is equal to this object, false otherwise.
  */
  bool equals(const Match * other) const {
    return (_tpos == other->_tpos &&
	    _bpos == other->_bpos &&
	    _l == other->_l &&
	    _n == other->_n &&
	    _orientation == other->_orientation);
  }


  /*!
  *   \fn const char * get_top_window(const char * top_forward) const
  *   \brief Returns a pointer into a character array. 
  *   \param top_forward A character array, representing the top sequence matched in this Match object (forward orientation).
  *   \return A pointer into the char array parameter top_forward at the top position from this Match.
  *   \sa _tpos
  *   
  *   Function simple returns a pointer into character array passed 
  *   as an argument, at the index indicated by _tpos (top position).
  */
  const char * get_top_window(const char * top_forward) const {
    return &top_forward[_tpos];
  }


  /*!
  *   \fn const char * get_bot_window(const char * bot_forward,
			      const char * bot_reverse,
			      const unsigned int bot_len) const
  *   \brief Returns a pointer into a character array.
  *   \param bot_forward A character array representing the bottom sequence of this Match with forward orientation.
  *   \param bot_reverse A character array representing the bottom sequence of this Match with reverse orientation.
  *   \param bot_len The length of the bottom sequence from this Match.
  *   \return A pointer into a character array representing the start of a Match in the bottom position.
  *   \sa _bpos
  *   
  *   Function takes two character arguments because the sequences are not stored locally in this Match object,
  *   only the positions.  The funtion wll return a pointer into the bot_reverse parameter if orientation is 
  *  equal to -1, otherwise it returns a pointer into bot_forward.
  */  
  const char * get_bot_window(const char * bot_forward,
			      const char * bot_reverse,
			      const unsigned int bot_len) const {
    if (_orientation == -1) {
      const unsigned int rev_pos = bot_len - (_bpos + 1);
      return &bot_reverse[rev_pos];
    } else {
      return &bot_forward[_bpos];
    }
  }
};

#endif // MATCH_HH
