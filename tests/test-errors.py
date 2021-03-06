#! /usr/bin/env python
import sys, os, gc

import paircomp, fasta

dirname = os.path.dirname(__file__)
thisdir = os.path.normpath(dirname) + '/'
bindir = os.path.join(thisdir, '../../bin/')

def test():
    el = fasta.load_single(thisdir + 'el.txt')
    br = fasta.load_single(thisdir + 'br.txt')
    re = fasta.load_single(thisdir + 're.txt')

    ### first break the algorithms.

    simple_cmp = None

    try:
        simple_cmp = paircomp.do_simple_nxn_compare("ATCCGRRRRUUUUU", br, 10, .5)
    except:
        pass
    assert simple_cmp is None

    # break on windowsizes
    try:
        simple_cmp = paircomp.do_simple_nxn_compare(el, br, -1, .5)
    except:
        pass
    assert simple_cmp is None

    try:
        simple_cmp = paircomp.do_simple_nxn_compare(el, br, 256, .5)
    except:
        pass
    assert simple_cmp is None

    try:
        # break on windowsizes in relation to sequence len
        simple_cmp = paircomp.do_simple_nxn_compare("AAAAAAAAAAAAAAAAAAA",
                                                    br, 20, .5)
    except:
        pass
    assert simple_cmp is None

    try:
        simple_cmp = paircomp.do_simple_nxn_compare(el, "AAAAAAAAAAAAAAAAAAA",
                                                    20, .5)
    except:
        pass
    assert simple_cmp is None

    # break on threshold
    try:
        simple_cmp = paircomp.do_simple_nxn_compare(el, br, 9, -.01)
    except:
        pass
    assert simple_cmp is None

    try:
        simple_cmp = paircomp.do_simple_nxn_compare(el, br, 9, 1.01)
    except:
        pass
    assert simple_cmp is None

    simple_cmp_1 = paircomp.do_simple_nxn_compare(el, br, 20, .5)

    simple_cmp_2 = paircomp.do_simple_nxn_compare(el, br, 22, .5)

    result = None
    try:
        result = simple_cmp_1.contains(simple_cmp_2)
    except:
        pass
    assert result is None

    try:
        result = simple_cmp_1.equals(simple_cmp_2)
    except:
        pass
    assert result is None

    try:
        result = simple_cmp_1.subtract(simple_cmp_2)
    except:
        pass
    assert result is None

    try:
        result = simple_cmp_1.intersect(simple_cmp_2)
    except:
        pass

    assert result is None

    ### test transitivity stuff

    windowsize = 20
    threshold = 0.7
    ab = paircomp.do_simple_nxn_compare(el, br, windowsize, threshold)
    ab2 = paircomp.do_simple_nxn_compare(el, br, windowsize+1, threshold)
    bc = paircomp.do_simple_nxn_compare(br, re, windowsize, threshold)
    bc2 = paircomp.do_simple_nxn_compare(br, re, windowsize+1, threshold)
    ac = paircomp.do_simple_nxn_compare(el, re, windowsize, threshold)

    results = None
    try:
        results = paircomp.build_transitive(ab, bc, el, br, threshold)
    except:
        pass
    assert results is None

    try:
        results = paircomp.build_transitive(ab, bc, br, re, threshold)
    except Exception, e:
        pass
    assert results is None

    try:
        results = paircomp.build_transitive(ab2, bc, el, re, threshold)
    except:
        pass
    assert results is None

    try:
        results = paircomp.build_transitive(ab, bc2, el, re, threshold)
    except:
        pass
    assert results is None

    try:
        results = paircomp.filter_transitively(ab2, bc, ac)
    except:
        pass
    assert results is None

    try:
        results = paircomp.filter_transitively(ab, bc2, ac)
    except:
        pass
    assert results is None

    print 'SUCCESS.'
