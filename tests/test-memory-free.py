#! /usr/bin/env python
import _testdir
import gc

import sys, os, fasta
import paircomp

assert paircomp.__version__ == "1.0"
print """
Testing version %s of the paircomp extension module, loaded from

\t%s
""" % (paircomp.__version__, _testdir.paircomp_build_dir,)

el = fasta.load_single('el.txt')
br = fasta.load_single('br.txt')
re = fasta.load_single('re.txt')

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
