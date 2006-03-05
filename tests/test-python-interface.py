#! /usr/bin/env python
import sys, os, gc

print '==> in %s' % (sys.argv[0],)
import _testdir
_testdir.paircomp_version_message()

import paircomp, fasta

def print_cmp(cmp):
    print '--'
    for i in range(0, cmp.top_len):
        l = cmp[i]
        if l:
            print '==>', i, l

el = fasta.load_single('el.txt')
br = fasta.load_single('br.txt')

### four different ways of doing analyses:

#
# execute paircomp
#

os.system('../bin/paircomp el.txt br.txt 20 .5 el-br-0.5.p.cmp >& /dev/null')
p_cmp = paircomp.parse_paircomp(open('el-br-0.5.p.cmp').read(),
                                len(el), len(br), 20)

#
# execute seqcomp
#

os.system('../bin/seqcomp el.txt br.txt 20 10 noxml el-br-0.5.s.cmp >& /dev/null')
s_cmp = paircomp.parse_seqcomp(open('el-br-0.5.s.cmp').read(),
                               len(el), len(br), 20)

#
# & use the python interface.
#

simple_cmp = paircomp.do_simple_nxn_compare(el, br, 20, .5)
rolling_cmp = paircomp.do_rolling_nxn_compare(el, br, 20, .5)

## also load in an old one

old_paircomp = paircomp.parse_paircomp(open('orig-el-br.txt').read(),
                                       len(el), len(br), 20)

## check to make sure they all return the same result!

print 'Testing equality of internal vs external analyses...',

assert p_cmp == s_cmp
assert p_cmp == simple_cmp
assert p_cmp == rolling_cmp
assert p_cmp == old_paircomp

print 'success'

# round-trip both formats:

print 'Testing save/load round-trip...',

simple_cmp.save('simple.cmp')
simple_cmp_2 = paircomp.parse_paircomp(open('simple.cmp').read(),
                                       len(el), len(br), 20)

assert simple_cmp == simple_cmp_2

simple_cmp.save_as_seqcomp('simple_seqcomp.cmp')
simple_cmp_3 = paircomp.parse_seqcomp(open('simple_seqcomp.cmp').read(),
                                      len(el), len(br), 20)
assert simple_cmp == simple_cmp_3

print 'success'

#
# ok, now abuse the Python interface.
#

print 'Testing basic functions...',

assert simple_cmp.windowsize == 20
assert simple_cmp.top_len == len(el)
assert simple_cmp.bot_len == len(br)

for i in range(0, len(el)):
    l = simple_cmp[i]
    if len(l):
        for m in l:
            assert m.top == i

print 'success'

print 'Testing reverse of top/bottom...',

a = simple_cmp.reverse_top()
b = a.reverse_top()
c = b.reverse_top()

assert simple_cmp == b
assert a == c
assert b != c

a = simple_cmp.reverse_bottom()
b = a.reverse_bottom()
c = b.reverse_bottom()

assert simple_cmp == b
assert a == c
assert b != c

el_rc = fasta.reverse_complement(el)
br_rc = fasta.reverse_complement(br)
simple_cmp_fr = paircomp.do_simple_nxn_compare(el, br_rc, 20, .5)
simple_cmp_rf = paircomp.do_simple_nxn_compare(el_rc, br, 20, .5)
simple_cmp_rr = paircomp.do_simple_nxn_compare(el_rc, br_rc, 20, .5)

assert simple_cmp == simple_cmp_fr.reverse_bottom()
assert simple_cmp == simple_cmp_rf.reverse_top()
assert simple_cmp == simple_cmp_rr.reverse_bottom().reverse_top()

print 'success'

print 'Testing inverse...',

a = simple_cmp.invert()
b = a.invert()
c = b.invert()

assert simple_cmp == b
assert a == c

simple_cmp_i = paircomp.do_simple_nxn_compare(br, el, 20, .5)
assert simple_cmp == simple_cmp_i.invert()

print 'success'

print 'Testing filter by threshold...',

# filter gradually...
a = simple_cmp.filter(.6)
b = simple_cmp.filter(.7)
c = simple_cmp.filter(.8)

# and all at once
d = simple_cmp.filter(.8)

assert d == c

print 'success'

###

print 'Testing set operations...',

a = b = simple_cmp
assert a.contains(b)
assert b.contains(a)

c = a - b
assert c.is_empty()
assert a.contains(c)
assert c == a - (a.intersect(b))

b = simple_cmp.filter(.6)
assert a.contains(b)
assert not b.contains(a)
c = a - b
assert not c.is_empty()

print 'success'

### 1bp conversion

print 'Testing 1bp conversion...',

a = b = simple_cmp
b = simple_cmp.filter(.6)
b = b.isolate_matching_bases(el, br)
a = a.isolate_matching_bases(el, br)

assert a.contains(b)

print 'success'

### test transitivity filter.

print 'Testing transitivity filtering...',

A='aaaaaaaaaaaaaaaaaaaa'
B='tttttttttttttttttttt'
C='tttttttttttttttttttt'

AB=paircomp.do_simple_nxn_compare(A, B, 20, .9)
BC=paircomp.do_simple_nxn_compare(B, C, 20, .9)
AC=paircomp.do_simple_nxn_compare(A, C, 20, .9)

l = [i.is_empty() for i in paircomp.parser.filter_transitively(AB, BC, AC) ]
assert l == [0, 0, 0]

A='tttttttttttttttttttt'
B='aaaaaaaaaaaaaaaaaaaa'
C='tttttttttttttttttttt'

AB=paircomp.do_simple_nxn_compare(A, B, 20, .9)
BC=paircomp.do_simple_nxn_compare(B, C, 20, .9)
AC=paircomp.do_simple_nxn_compare(A, C, 20, .9)

l = [i.is_empty() for i in paircomp.parser.filter_transitively(AB, BC, AC) ]
assert l == [0, 0, 0]

A='tttttttttttttttttttt'
B='tttttttttttttttttttt'
C='aaaaaaaaaaaaaaaaaaaa'

AB=paircomp.do_simple_nxn_compare(A, B, 20, .9)
BC=paircomp.do_simple_nxn_compare(B, C, 20, .9)
AC=paircomp.do_simple_nxn_compare(A, C, 20, .9)

l = [i.is_empty() for i in paircomp.parser.filter_transitively(AB, BC, AC) ]
assert l == [0, 0, 0]

### try some broken ones.

A='ttttttttttttttttttttt'
B='tttttttttttttttttttt'
C='aaaaaaaaaaaaaaaaaaa'

AB=paircomp.do_simple_nxn_compare(A, B, 15, .9)
BC=paircomp.do_simple_nxn_compare(B, C, 15, .9)
AC=paircomp.do_simple_nxn_compare(A, C, 15, .9)

failed = False
try:
    paircomp.parser.filter_transitively(AB, BC, AB) # should fail: AB-->AC
    #    ---------------------------------------^^
except:
    failed = True

assert failed

#AB=paircomp.do_simple_nxn_compare(A, B, 20, .9)
#BC=paircomp.do_simple_nxn_compare(B, C, 20, .9)
#AC=paircomp.do_simple_nxn_compare(A, C, 19, .9)
#    -----------------------------------^^

#failed = False
#try:
#    paircomp.parser.filter_transitively(AB, BC, AC) # should fail: wrong wsize
#except:
#    failed = True

#assert failed

print 'success'

###

print '\nGeneral tests COMPLETE.  Enjoy your day!'

gc.collect()                            # for valgrind-style stuff.
