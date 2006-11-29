from py_nway import NwayComparison as PyNwayComparison, do_3way_regression
from paircomp.nway import NwayComparison

def check_result(seqs, result):
    py_nway_cmp = PyNwayComparison(10, 1.0, *seqs)
    nway_cmp = NwayComparison(10, 1.0, *seqs)

    p = py_nway_cmp.filter()
    assert repr(p) == result

    p2 = nway_cmp.get_filtered_paths()
    p2.sort()

    assert p == p2

    if len(seqs) == 3:
        p3 = do_3way_regression(seqs[0], seqs[1], seqs[2], 10, 1.0)
        p3.sort()
        
        assert p2 == p3

def test():

    seq = 'AAAAAAAAAA'
    seq2 = 'TTTTTTTTTT'

    trials = (

        ( (seq, seq, seq), '[[0(+), 0(+), 0(+)]]' ),
        ( (seq, seq, seq2), '[[0(+), 0(-), 0(+)]]' ),
        ( (seq, seq2, seq), '[[0(-), 0(-), 0(+)]]' ),
        ( (seq, seq2, seq2), '[[0(-), 0(+), 0(+)]]' ),

        ( (seq, seq, seq, seq), '[[0(+), 0(+), 0(+), 0(+)]]' ),
        ( (seq, seq, seq, seq2), '[[0(+), 0(+), 0(-), 0(+)]]' ),
        ( (seq, seq, seq2, seq2), '[[0(+), 0(-), 0(+), 0(+)]]' ),
        ( (seq, seq2, seq2, seq2), '[[0(-), 0(+), 0(+), 0(+)]]' ),
        ( (seq, seq2, seq, seq), '[[0(-), 0(-), 0(+), 0(+)]]' ),
        ( (seq, seq, seq2, seq), '[[0(+), 0(-), 0(-), 0(+)]]' ),
        ( (seq, seq2, seq, seq2), '[[0(-), 0(-), 0(-), 0(+)]]' ),
        ( (seq, seq2, seq2, seq), '[[0(-), 0(+), 0(-), 0(+)]]' ),

        ( (seq, seq2, seq, seq, seq2), '[[0(-), 0(-), 0(+), 0(-), 0(+)]]' ),

        )

    for (seqs, result) in trials:
        check_result(seqs, result)

def test2():
    
    seq = 'AAAAAAAAAAA'             # 11

    nway_cmp = NwayComparison(10, 1.0, seq, seq, seq)

    p = nway_cmp.get_filtered_paths()
    assert repr(p) == '[[0(+), 0(+), 0(+)], [0(+), 0(+), 1(+)], [0(+), 1(+), 1(+)], [0(+), 1(+), 0(+)], [1(+), 1(+), 1(+)], [1(+), 1(+), 0(+)], [1(+), 0(+), 0(+)], [1(+), 0(+), 1(+)]]'

    p2 = do_3way_regression(seq, seq, seq, 10, 1.0)

    p.sort()
    p2.sort()

    assert p == p2

def test3():

    seq = 'AAAAAAAAAAA'             # 11
    seq2 = 'TTTTTTTTTTT'            # 11

    nway_cmp = NwayComparison(10, 1.0, seq, seq2, seq)

    p = nway_cmp.get_filtered_paths()
    assert repr(p) == '[[0(-), 1(-), 0(+)], [0(-), 1(-), 1(+)], [0(-), 0(-), 1(+)], [0(-), 0(-), 0(+)], [1(-), 0(-), 1(+)], [1(-), 0(-), 0(+)], [1(-), 1(-), 0(+)], [1(-), 1(-), 1(+)]]'

    p2 = do_3way_regression(seq, seq2, seq, 10, 1.0)

    p.sort()
    p2.sort()

    assert p == p2
