#! /usr/bin/env python
import sys
sys.path.append('../python/')
import os

import fasta                            # from cartwheel
from paircomp import do_rolling_nxn_compare

# construct a quick seqcomp
el = fasta.load_single('el.txt')
br = fasta.load_single('br.txt')

a = do_rolling_nxn_compare(el, br, 50, .5)

b = a.isolate_matching_bases(el, br)

print 'test-1bp-matches.py PASSED.'
