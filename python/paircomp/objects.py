import _paircomp_parser
import string

#
# Comparison
#

class Comparison:
    """
    paircomp Comparison object.

    Attributes:
      * windowsize, top_len, bot_len

    Methods:
      * __init__(cmp, top_len, bot_len, windowsize)
      * save(filename)
      * save_as_seqcomp(filename)
      * reverse_top()
      * reverse_bottom()
      * filter(float_threshold)
      * contains(other_cmp)
      * equals(other_cmp), a.k.a '==' operator
      * subtract(other_cmp), a.k.a. '-' operator
      * is_empty()
      * get_matches, a.k.a. []
    """
    
    def __init__(self, _c_object):
        """
        'source' is a string containing a paircomp-format comparison.
        """
        self._comparisonObj = _c_object
        self.windowsize = _paircomp_parser.get_windowsize(_c_object)
        self.top_len = _paircomp_parser.get_top_length(_c_object)
        self.bot_len = _paircomp_parser.get_bottom_length(_c_object)

    def save(self, filename):
        _paircomp_parser.save_paircomp_comparison(filename, self._comparisonObj)

    def save_as_seqcomp(self, filename):
        _paircomp_parser.save_seqcomp_comparison(filename, self._comparisonObj)

    def reverse_top(self):
        r = _paircomp_parser.reverse_top(self._comparisonObj)
        return Comparison(r)

    def reverse_bottom(self):
        r = _paircomp_parser.reverse_bot(self._comparisonObj)
        return Comparison(r)

    def invert(self):
        r = _paircomp_parser.invert(self._comparisonObj)
        return Comparison(r)

    def intersect(self, other):
        assert isinstance(other, Comparison)
        r = _paircomp_parser.intersect(self._comparisonObj, other._comparisonObj)
        return Comparison(r)

    def filter(self, threshold):
        f = _paircomp_parser.filter_matches(self._comparisonObj, threshold)
        return Comparison(f)

    def isolate_matching_bases(self, top_seq, bot_seq):
        f = _paircomp_parser.isolate_matching_bases(self._comparisonObj,
                                                    top_seq, bot_seq)
        return Comparison(f)

    def contains(self, other):
        assert isinstance(other, Comparison)
        b = _paircomp_parser.contains(self._comparisonObj,
                                      other._comparisonObj)
        return b

    def equals(self, other):
        if self.contains(other) and other.contains(self):
            return 1

        return 0

    def __eq__(self, other):
        return self.equals(other)

    def __ne__(self, other):
        return not self.equals(other)

    def subtract(self, other):
        cmp = _paircomp_parser.subtract(self._comparisonObj,
                                        other._comparisonObj)
        return Comparison(cmp)

    def __sub__(self, other):
        return self.subtract(other)

    def is_empty(self):
        return _paircomp_parser.is_empty(self._comparisonObj)

    def get_matches(self, i):
        l = _paircomp_parser.get_matches(self._comparisonObj, i)

        return [ Match(a, b, c, d, e) for (a, b, c, d, e) in l ]
        
        return MatchList(l)

    def __getitem__(self, i):
        return self.get_matches(i)

#
# Match
#

class Match:
    """
    An individual match.  Contains top, bot, length, matches, and orientation.
    """
    def __init__(self, top, bot, length, matches, orientation):
        self.top = top
        self.bot = bot
        self.length = length
        self.matches = matches
        self.orientation = orientation

    def __repr__(self):
        return "[ %d --> %d, %d/%d matches (orientation %d) ]" % (self.top,
                                                              self.bot,
                                                              self.matches,
                                                              self.length,
                                                              self.orientation)
