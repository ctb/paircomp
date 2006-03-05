#! /usr/bin/env python
"""
A very quick & dirty sister-output-maker for seqcomp-style comparisons.
"""
import sys
from biolib import fasta

if len(sys.argv) != 6:
    sys.stderr.write("""\
Usage:

   %s seqfile1 seqfile2 windowsize comparison outfile

""" % (sys.argv[0],))
    sys.exit(-1)

windowsize = int(sys.argv[3])

#
# Open everything up & read it.
#

seq1_d = fasta.load(open(sys.argv[1]))
assert len(seq1_d.keys()) == 1, "can only have one sequence in %s!" % (sys.argv[1],)

seq1_name, seq1_val = seq1_d.items()[0]

seq2_d = fasta.load(open(sys.argv[2]))
assert len(seq2_d.keys()) == 1, "can only have one sequence in %s!" % (sys.argv[2],)

seq2_name, seq2_val = seq2_d.items()[0]

cmp = open(sys.argv[4]).read()
outfile = open(sys.argv[5], 'w')

# here is where we would do minimal consistency checking, if we were
# brave.  we're not.

outfile.write("""\
<Start>
<First>
  <Sequence name="%s">
%s
  </Sequence>
</First>
<Second>
  <Sequence name="%s">
%s
  </Sequence>
</Second>
<Comparison type="shoudan" name="cmp" topStart="0" botStart="0" topEnd="%d" botEnd="%d"
 windowsize="%d" minThreshold="0.0">
%s
</Comparison>
</Start>
""" % (seq1_name, seq1_val,
       seq2_name, seq2_val,
       len(seq1_val),
       len(seq2_val),
       windowsize,
       cmp))
