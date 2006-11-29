#! /usr/bin/env python
import sys, os, gc

print '==> in %s' % (sys.argv[0],)
import _testdir
_testdir.paircomp_version_message()

import paircomp, fasta

A = "A"*21
B = "T"*21

el = fasta.load_single('el.txt')
br = fasta.load_single('br.txt')
re = fasta.load_single('re.txt')

def print_cmp(cmp):
    print '---'
    for i in range(0, cmp.top_len):
        l = cmp[i]
        if l:
            print '==>', i, l

def test_transitive(A, B, C, windowsize=20, threshold=.5):
    ab = paircomp.do_simple_nxn_compare(A, B, windowsize, threshold)
    bc = paircomp.do_simple_nxn_compare(B, C, windowsize, threshold)
    ac = paircomp.do_simple_nxn_compare(A, C, windowsize, threshold)

    ac_tran1 = paircomp.build_transitive(ab, bc, A, C, threshold)
    ac_tran2 = paircomp.filter_transitively(ab, bc, ac)[2]

    assert ac_tran1 == ac_tran2

test_transitive(A, A, A)
test_transitive(A, A, B)

test_transitive(A, B, A)
test_transitive(A, B, B)

test_transitive(B, A, A)
test_transitive(B, A, B)

test_transitive(B, B, A)
test_transitive(B, B, B)

print 'el-br-re'

test_transitive(el, br, re)
test_transitive(el, fasta.rc(br), re)
test_transitive(el, br, fasta.rc(re))
test_transitive(el, fasta.rc(br), fasta.rc(re))

test_transitive(fasta.rc(el), br, re)
test_transitive(fasta.rc(el), fasta.rc(br), re)
test_transitive(fasta.rc(el), br, fasta.rc(re))
test_transitive(fasta.rc(el), fasta.rc(br), fasta.rc(re))

sys.exit(0)
