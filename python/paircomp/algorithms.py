import _paircomp_algorithms
from objects import Comparison

def do_simple_nxn_compare(seq1, seq2, windowsize, threshold):

    cmp = _paircomp_algorithms.do_simple_nxn_compare(seq1, seq2, windowsize,
                                                     threshold)

    return Comparison(cmp)


def do_rolling_nxn_compare(seq1, seq2, windowsize, threshold):

    cmp = _paircomp_algorithms.do_rolling_nxn_compare(seq1, seq2, windowsize,
                                                      threshold)

    return Comparison(cmp)

def do_hashed_n_compare(seq1, seq2, windowsize, threshold):

    cmp = _paircomp_algorithms.do_hashed_n_compare(seq1, seq2, windowsize,
                                                   threshold)

    return Comparison(cmp)
