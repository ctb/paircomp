#! /usr/bin/env python
import sys
from paircomp import parse_seqcomp

cmp_txt = open(sys.argv[1]).read()
cmp = parse_seqcomp(cmp_txt)
cmp.save(sys.argv[2])
