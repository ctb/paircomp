#! /usr/bin/env python
import sys

import paircomp

def test():
    A='aaaaaaaaaaaaaaaaaaaa'
    B='aaaaaaaaaaaaaaaaaaaa'
    C='aaaaaaaaaaaaaaaaaaaa'

    AB=paircomp.do_simple_nxn_compare(A, B, 20, .9)
    BC=paircomp.do_simple_nxn_compare(B, C, 20, .9)
    AC=paircomp.do_simple_nxn_compare(A, C, 20, .9)

    l = [i.is_empty() for i in paircomp.parser.filter_transitively(AB, BC, AC) ]

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

    print 'Tests PASSED.'
