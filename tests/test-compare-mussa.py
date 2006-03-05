#! /usr/bin/env python
import sys, os, gc

print '==> in %s' % (sys.argv[0],)
import _testdir
_testdir.paircomp_version_message()

print _testdir.paircomp_build_dir

import paircomp, fasta

def print_cmp(cmp):
    print '--'
    for i in range(0, cmp.top_len):
        l = cmp[i]
        if l:
            print '==>', i, l

t1 = fasta.load_single('t1.fa', force='DNA')
t2 = fasta.load_single('t2.fa', force='DNA')
t3 = fasta.load_single('t3.fa', force='DNA')

AB = paircomp.do_simple_nxn_compare(t1, t2, 4, 1.0)
BC = paircomp.do_simple_nxn_compare(t2, t3, 4, 1.0)
AC = paircomp.do_simple_nxn_compare(t1, t3, 4, 1.0)

(new_AB, new_BC, new_AC) = paircomp.filter_transitively(AB, BC, AC)

#print_cmp(new_AB)
#print_cmp(new_BC)
#print_cmp(new_AC)

expected = \
"""
4:0,0,20
4:0,0,-5
4:0,29,20
4:0,29,-5
4:0,-16,-5
4:0,-16,20
4:33,6,33
4:33,6,-33
4:33,43,33
4:33,43,-33
4:33,-43,-33
4:33,-43,33
4:33,-6,-33
4:33,-6,33
"""

def parse_tristan_3way(s):
    """
    Parse Tristan's output.
    """
    l = []
    for line in s.split('\n'):
        if line:
            w, rest = line.split(':')
            a, b, c = rest.split(',')
            w, a, b, c = map(int, (w, a, b, c,))

            l.append((a, b, c,))

    return l

def print_tristan_3way(ab, bc, ac):
    """
    Convert my ab/bc/ac transitive into Tristan's format.
    """
    w = ab.windowsize
    l = []
    for i in range(0, ab.top_len):
        top_matches = ab[i]
        for m1 in top_matches:
            c_matches = ac[m1.top]
            for m2 in c_matches:

                # reverse middle?
                if m1.orientation == -1:
                    b = -(m1.bot - (w - 1))
                else:
                    b = m1.bot

                # reverse bottom?
                if m2.orientation == -1:
                    c = -(m2.bot - (w - 1))
                else:
                    c = m2.bot
                    
                l.append((m1.top, b, c,))

    return l

print 'Comparing paircomp output to Mussa\'s...',

l = print_tristan_3way(new_AB, new_BC, new_AC)
l.sort()

l2 = parse_tristan_3way(expected)
l2.sort()

assert l == l2
print 'success'

if l != l2:
    print 'DIFF: CTB'

    for x in l:
        if not x in l2:
            print x

    print 'vs TDB'

    for x in l2:
        if not x in l:
            print x
