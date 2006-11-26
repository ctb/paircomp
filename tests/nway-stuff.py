#! /usr/bin/env python
import sys
from biolib import fasta

####

import _testdir
import nway
import paircomp

seqs = [ fasta.load_single(open(f)) for f in sys.argv[1:] ]

print 'loaded %d sequences' % (len(seqs),)

nway_cmp = nway.NwayComparison(20, 0.7, seqs[0])
for seq in seqs[1:]:
    nway_cmp.add_sequence(seq)
