# See README.txt for license and copyright information.

"""
File parsers and objects for seqcomp and shoudan-style sequence comparisons.

(This module is a Python wrapper around the _paircomp_parser.so C/C++ extension
module.)
"""

import _paircomp_parser
import string
from objects import Comparison

def parse_paircomp(buf, top_length, bot_length, windowsize):
    _c_obj = _paircomp_parser.parse_paircomp_comparison(buf,
                                                        top_length, bot_length,
                                                        windowsize)
    return Comparison(_c_obj)

def parse_seqcomp(buf, top_length, bot_length, windowsize):
    _c_obj = _paircomp_parser.parse_seqcomp_comparison(buf,
                                                       top_length,
                                                       bot_length,
                                                       windowsize)
    
    return Comparison(_c_obj)

def build_transitive(ab, bc, seq_a, seq_c, threshold):
    _c_obj = _paircomp_parser.build_transitive(ab._comparisonObj,
                                               bc._comparisonObj,
                                               seq_a, seq_c, threshold)
    return Comparison(_c_obj)

def filter_transitively(ab, bc, ac):
    
    l = _paircomp_parser.filter_transitively(ab._comparisonObj,
                                             bc._comparisonObj,
                                             ac._comparisonObj)
    (new_ab, new_bc, new_ac) = l

    return (Comparison(new_ab),
            Comparison(new_bc),
            Comparison(new_ac))
