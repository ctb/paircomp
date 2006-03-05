#! /usr/bin/env python
import sys, os, fasta

print '==> in %s' % (sys.argv[0],)
import _testdir
_testdir.paircomp_version_message()

import paircomp, fasta

files = ['delta', 'gcm', 'otx']
types = ['v0.8', 'v0.9']

for file in files:
    seq1 = fasta.load_single('regression-data/%s-sp.txt' % (file,))
    seq2 = fasta.load_single('regression-data/%s-lv.txt' % (file,))

    ext1 = types[0]
    filename1 = 'regression-data/%s-%s.cmp' % (file, ext1,)
    
    cmp1 = open(filename1).read()
    cmp1 = paircomp.parse_paircomp(cmp1, len(seq1), len(seq2), 20)

    ext2 = types[1]
    filename2 = 'regression-data/%s-%s.cmp' % (file, ext1,)

    cmp2 = open(filename2).read()
    cmp2 = paircomp.parse_paircomp(cmp2, len(seq1), len(seq2), 20)

    print 'Comparing %s-%s to %s-%s...' % (file, ext1, file, ext2,),

    assert(cmp1 == cmp2)
    print 'match!'

print '\nRegression tests COMPLETE.  Enjoy your day!'
