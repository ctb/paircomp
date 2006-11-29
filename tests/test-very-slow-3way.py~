#! /usr/bin/env python
import sys

print '==> in %s' % (sys.argv[0],)
import _testdir
_testdir.paircomp_version_message()

import paircomp
import time
import sys
import fasta

thresh = 14.0/21

A='tttttttttttttttttttttttttttttttttttttttttttttttttt'
B='tttttttttttttttttttttttttttttttttttttttttttttttttt'
C='tttttttttttttttttttttttttttttttttttttttttttttttttt'

AB=paircomp.do_rolling_nxn_compare(A, B, 21, thresh)
BC=paircomp.do_rolling_nxn_compare(B, C, 21, thresh)
AC=paircomp.do_rolling_nxn_compare(A, C, 21, thresh)

a=time.time()
tran=paircomp.parser.filter_transitively(AB,BC,AC)
b=time.time()
print 'A B C are same '+str(b-a)
