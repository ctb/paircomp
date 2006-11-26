#! /usr/bin/env python
import sys
from biolib import fasta

####

import _testdir
import nway
import paircomp

seqs = [ fasta.load_single(open(f)) for f in sys.argv[1:] ]

print 'loaded %d sequences' % (len(seqs),)

nway_cmp = nway.NwayComparison(20, 0.7, *seqs)

p = nway_cmp.filter()
p.sort()

if len(seqs) == 3:
    p2 = nway.do_3way_regression(seqs[0], seqs[1], seqs[2], 20, 0.7)
    p2.sort()

#assert p == p2
for i in p:
    if i not in p2:
        print i

for i in p2:
    if i not in p:
        print i
        (x, y, z) = i
        
        a = seqs[0][x.pos:x.pos+20]
        b = seqs[1][y.pos:y.pos+20]
        c = seqs[2][z.pos:z.pos+20]

        Nab = 0
        for (m, n) in zip(a, b):
            if m == n: Nab += 1

        Nbc = 0
        for (m, n) in zip(b, c):
            if m == n: Nbc += 1

        Nac = 0
        for (m, n) in zip(a, c):
            if m == n: Nac += 1

        print Nab, Nbc, Nac
