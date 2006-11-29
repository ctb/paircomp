#! /usr/bin/env python
import os
import fasta

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname) + '/'
bindir = os.path.join(thisdir, '../../bin/')

from paircomp import do_rolling_nxn_compare

def test():
    # construct a quick seqcomp
    el = fasta.load_single(thisdir + 'el.txt')
    br = fasta.load_single(thisdir + 'br.txt')

    a = do_rolling_nxn_compare(el, br, 50, .5)

    b = a.isolate_matching_bases(el, br)

    print 'test-1bp-matches.py PASSED.'
