#! /usr/bin/env python
import sys

## import local Fasta-format library.
sys.path.append('../tests/')
import fasta

if len(sys.argv) != 6:
    sys.stderr.write("""Usage:

    %s top_seq bot_seq paircomp_file windowsize threshold
    
""" % (sys.argv[0],))

    sys.exit(-1)

top_seq = fasta.load_single(sys.argv[1])
bot_seq = fasta.load_single(sys.argv[2])
paircomp_file = open(sys.argv[3])
windowsize = int(sys.argv[4])
threshold = float(sys.argv[5])

print """
<Start>
 <First>
  <Sequence name="top sequence">
%s
  </Sequence>
  
  <!-- top sequence annotations can go here -->
  
 </First>
 <Second>
  <Sequence name="bottom sequence">
%s
  </Sequence>
  
  <!-- top sequence annotations can go here -->
  
 </Second>

 <!-- One or more comparisons -->

 <Comparison name="a comparison" type="shoudan" topStart="0" botStart="0"
        topLength="%d" windowsize="%d" minThreshold="%f">
%s
 </Comparison>
</Start>
""" % (top_seq, bot_seq, len(top_seq) - windowsize + 1, windowsize, threshold,
       paircomp_file.read())
