#! /usr/bin/env python
import gc

import sys, os, fasta
import paircomp

assert paircomp.__version__ == "1.0"

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname) + '/'
bindir = os.path.join(thisdir, '../../bin/')

def test():

    el = fasta.load_single(thisdir + 'el.txt')
    br = fasta.load_single(thisdir + 'br.txt')
    re = fasta.load_single(thisdir + 're.txt')

    print 'testing comparison:'
    print 'GC\tloop count'

    for i in range(0, 100):
        if i % 10 == 0:
            print gc.collect(), '\t', i
        cmp = paircomp.do_rolling_nxn_compare(el, br, 20, .5)

    ### test transitivity, too.

    ab = paircomp.do_rolling_nxn_compare(el, br, 20, .7)
    bc = paircomp.do_rolling_nxn_compare(br, re, 20, .7)
    ac = paircomp.do_rolling_nxn_compare(el, re, 20, .7)

    print 'testing transitivity:'
    print 'GC\tloop count'
    for i in range(0, 100):
        if i % 10 == 0:
            print gc.collect(), '\t', i
        (new_ab, new_bc, new_ac) = paircomp.filter_transitively(ab, bc, ac)

    del new_ab, new_bc, new_ac
    del ab, bc, ac

    print 'remaining objects to collect:'
    print gc.collect()

    # there should be no extraneous memory left over, except for stuff used
    # by Python itself.
