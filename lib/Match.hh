// Class definition for fixed-width window match.
//
// See README.txt for license and copyright information.

#ifndef MATCH_HH		// guard against multiple inclusion
#define MATCH_HH

#include <assert.h>

class Match {
protected:
  unsigned int _tpos, _bpos;
  unsigned char _l, _n;
  char _orientation;
public:
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

  unsigned int get_top_pos() const { return _tpos; };
  unsigned int get_bot_pos() const { return _bpos; };
  unsigned int get_length() const { return (unsigned int) _l; };
  unsigned int get_n_matching() const { return (unsigned int) _n; };
  int get_orientation() const { return (int) _orientation; };

  Match * copy() const {
    return new Match(_tpos, _bpos, _l, _n, _orientation);
  }

  void reverse_top(unsigned int top_len, unsigned int bot_len) {
    if (_orientation == -1) {	// good gawd
      _bpos = _bpos - (_l - 1);
    } else {
      _bpos = _bpos + (_l - 1);
    }

    _tpos = top_len - _tpos - _l;
    _orientation = -1 * _orientation;
  }

  void reverse_bottom(unsigned int bot_len) {
    _bpos = bot_len - (_bpos + 1);
    _orientation = -1 * _orientation;
  }

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

  bool equals(const Match * other) const {
    return (_tpos == other->_tpos &&
	    _bpos == other->_bpos &&
	    _l == other->_l &&
	    _n == other->_n &&
	    _orientation == other->_orientation);
  }

  const char * get_top_window(const char * top_forward) const {
    return &top_forward[_tpos];
  }
  
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
