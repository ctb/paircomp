#! /usr/bin/env python
import sys, os, gc

import paircomp, fasta

A = "A"*21
B = "T"*21

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname) + '/'
bindir = os.path.join(thisdir, '../../bin/')

el = fasta.load_single(thisdir + 'el.txt')
br = fasta.load_single(thisdir + 'br.txt')
re = fasta.load_single(thisdir + 're.txt')

def is_transitive(A, B, C, windowsize=20, threshold=.5):
    ab = paircomp.do_simple_nxn_compare(A, B, windowsize, threshold)
    bc = paircomp.do_simple_nxn_compare(B, C, windowsize, threshold)
    ac = paircomp.do_simple_nxn_compare(A, C, windowsize, threshold)

    ac_tran1 = paircomp.build_transitive(ab, bc, A, C, threshold)
    ac_tran2 = paircomp.filter_transitively(ab, bc, ac)[2]

    assert ac_tran1 == ac_tran2

def test():
    is_transitive(A, A, A)
    is_transitive(A, A, B)

    is_transitive(A, B, A)
    is_transitive(A, B, B)

    is_transitive(B, A, A)
    is_transitive(B, A, B)

    is_transitive(B, B, A)
    is_transitive(B, B, B)

    print 'el-br-re'

    is_transitive(el, br, re)
    is_transitive(el, fasta.rc(br), re)
    is_transitive(el, br, fasta.rc(re))
    is_transitive(el, fasta.rc(br), fasta.rc(re))

    is_transitive(fasta.rc(el), br, re)
    is_transitive(fasta.rc(el), fasta.rc(br), re)
    is_transitive(fasta.rc(el), br, fasta.rc(re))
    is_transitive(fasta.rc(el), fasta.rc(br), fasta.rc(re))
