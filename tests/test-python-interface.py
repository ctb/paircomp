#! /usr/bin/env python
import sys, os, gc

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname)
bindir = os.path.join(thisdir, '../bin/')

import paircomp, fasta

class Test:
    def setup(self):
        self.el_path = os.path.join(thisdir, 'el.txt')
        self.br_path = os.path.join(thisdir, 'br.txt')

        self.el = fasta.load_single(self.el_path)
        self.br = fasta.load_single(self.br_path)

        self.cmp = paircomp.do_rolling_nxn_compare(self.el, self.br, 20, .5)

    def test(self):
        el = self.el
        br = self.br

        el_path = self.el_path
        br_path = self.br_path

        ### four different ways of doing analyses:

        #
        # execute paircomp
        #

        outname = 'el-br-0.5.p.cmp'
        os.system('rm %s >& /dev/null' % (outname,))
        os.system('%s/paircomp %s %s 20 .5 %s >& /dev/null' % \
                  (bindir, self.el_path, self.br_path, outname))
    
        p_cmp = paircomp.parse_paircomp(open(outname).read(),
                                        len(el), len(br), 20)
        os.system('rm %s >& /dev/null' % (outname,))
        
        #
        # execute seqcomp
        #

        outname = 'el-br-0.5.s.cmp'
        os.system('rm %s >& /dev/null' % (outname,))
        os.system('%s/seqcomp %s %s 20 10 noxml %s >& /dev/null' % \
                  (bindir, el_path, br_path, outname))
        s_cmp = paircomp.parse_seqcomp(open(outname).read(),
                                       len(el), len(br), 20)
        os.system('rm %s >& /dev/null' % (outname,))

        #
        # & use the python interface.
        #

        rolling_cmp = self.cmp
        simple_cmp = paircomp.do_simple_nxn_compare(el, br, 20, .5)

        ## also load in an old analysis for regression testing

        old_paircomp = paircomp.parse_paircomp(open('%s/orig-el-br.txt' % \
                                                    (thisdir,)).read(),
                                               len(el), len(br), 20)

        ## check to make sure they all return the same result!

        print 'Testing equality of internal vs external analyses...',

        assert p_cmp == s_cmp
        assert p_cmp == simple_cmp
        assert p_cmp == rolling_cmp
        assert p_cmp == old_paircomp

    def test_round_trip(self):
        "round-trip both formats"

        el, br = self.el, self.br
        cmp = self.cmp

        print 'Testing save/load round-trip...'

        cmp.save('simple.cmp')
        simple_cmp_2 = paircomp.parse_paircomp(open('simple.cmp').read(),
                                               len(el), len(br), 20)

        assert cmp == simple_cmp_2
        os.system('rm simple.cmp')

        cmp.save_as_seqcomp('simple_seqcomp.cmp')
        r = open('simple_seqcomp.cmp').read()
        simple_cmp_3 = paircomp.parse_seqcomp(r, len(el), len(br), 20)
        assert cmp == simple_cmp_3

    def test_python_functions(self):

        el, br = self.el, self.br
        cmp = self.cmp
        
        #
        # ok, now abuse the Python interface.
        #

        print 'Testing basic functions...',

        assert cmp.windowsize == 20
        assert cmp.top_len == len(el)
        assert cmp.bot_len == len(br)

        n = []
        for i in range(0, len(el)):
            l = cmp[i]
            for m in l:
                assert m.top == i
            n.extend(l)

        n2 = []
        for match in cmp:
            n2.append(match)

        for (a, b) in zip(n, n2):
            assert a == b

        print 'Testing reverse of top/bottom...',

        a = cmp.reverse_top()
        b = a.reverse_top()
        c = b.reverse_top()

        assert cmp == b
        assert a == c
        assert b != c

        a = cmp.reverse_bottom()
        b = a.reverse_bottom()
        c = b.reverse_bottom()

        assert cmp == b
        assert a == c
        assert b != c

        el_rc = fasta.reverse_complement(el)
        br_rc = fasta.reverse_complement(br)
        cmp_fr = paircomp.do_simple_nxn_compare(el, br_rc, 20, .5)
        cmp_rf = paircomp.do_simple_nxn_compare(el_rc, br, 20, .5)
        cmp_rr = paircomp.do_simple_nxn_compare(el_rc, br_rc, 20, .5)
        
        assert cmp == cmp_fr.reverse_bottom()
        assert cmp == cmp_rf.reverse_top()
        assert cmp == cmp_rr.reverse_bottom().reverse_top()

        print 'Testing inverse...',

        a = cmp.invert()
        b = a.invert()
        c = b.invert()

        assert cmp == b
        assert a == c

        cmp_i = paircomp.do_simple_nxn_compare(br, el, 20, .5)
        assert cmp == cmp_i.invert()
            
        print 'Testing filter by orientation...',

        # should be equal
        a = cmp.filter_orientation(True, True)
        assert a == cmp

        # should not be equal!
        a = cmp.filter_orientation(True, False)
        assert a != cmp

        b = cmp.filter_orientation(False, True)
        assert a != cmp
        
        # should be equal
        a = cmp.filter_orientation(True, False)
        b = cmp.reverse_bottom().filter_orientation(False, True)
        b = b.reverse_bottom()

        assert a == b

        # should be equal
        a = cmp.filter_orientation(True, False)
        b = cmp.reverse_top().filter_orientation(False, True)
        b = b.reverse_top()

        assert a == b

        print 'Testing filter by threshold...',
            
        # filter gradually...
        a = cmp.filter(.6)
        b = cmp.filter(.7)
        c = cmp.filter(.8)

        # and all at once
        d = cmp.filter(.8)

        assert d == c

        ###

        print 'Testing set operations...',

        a = b = cmp
        assert a.contains(b)
        assert b.contains(a)

        c = a - b
        assert c.is_empty()
        assert a.contains(c)
        assert c == a - (a.intersect(b))

        b = cmp.filter(.6)
        assert a.contains(b)
        assert not b.contains(a)
        c = a - b
        assert not c.is_empty()

        ### 1bp conversion

        print 'Testing 1bp conversion...',

        a = b = cmp
        b = cmp.filter(.6)
        b = b.isolate_matching_bases(el, br)
        a = a.isolate_matching_bases(el, br)

        assert a.contains(b)

    def test_3way(self):
        "Test transitivity filtering"

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
            # should fail: AB-->AC
            paircomp.parser.filter_transitively(AB, BC, AB)
            #    ---------------------------------------^^
        except:
            failed = True

        assert failed
