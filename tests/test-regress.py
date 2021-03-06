#! /usr/bin/env python
import sys, os, fasta

import paircomp, fasta

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname) + '/'
bindir = os.path.join(thisdir, '../../bin/')

files = ['delta', 'gcm', 'otx']
types = ['v0.8', 'v0.9']

def test():
    for file in files:
        seq1 = fasta.load_single(thisdir + 'regression-data/%s-sp.txt' % (file,))
        seq2 = fasta.load_single(thisdir + 'regression-data/%s-lv.txt' % (file,))

        ext1 = types[0]
        filename1 = thisdir + 'regression-data/%s-%s.cmp' % (file, ext1,)

        cmp1 = open(filename1).read()
        cmp1 = paircomp.parse_paircomp(cmp1, len(seq1), len(seq2), 20)

        ext2 = types[1]
        filename2 = thisdir + 'regression-data/%s-%s.cmp' % (file, ext1,)

        cmp2 = open(filename2).read()
        cmp2 = paircomp.parse_paircomp(cmp2, len(seq1), len(seq2), 20)

        print 'Comparing %s-%s to %s-%s...' % (file, ext1, file, ext2,),

        assert(cmp1 == cmp2)
        print 'match!'

    print '\nRegression tests COMPLETE.  Enjoy your day!'
