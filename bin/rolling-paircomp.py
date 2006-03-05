#! /usr/bin/env python
import sys
import fasta
import paircomp

top = fasta.load_single(sys.argv[1])
bot = fasta.load_single(sys.argv[2])
windowsize = int(sys.argv[3])
threshold = float(sys.argv[4])
outfile = sys.argv[5]

cmp = paircomp.do_rolling_nxn_compare(top, bot, windowsize, threshold)
cmp.save(outfile)
