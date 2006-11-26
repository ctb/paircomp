import _testdir
import paircomp

class NwayComparison:
    """
    Python implementation of the N-way transitivity algorithm used by Mussa.

    I've written an implementation that I think I understand, but no
    guarantees!

    Consider a four-way comparison.  The following paircomps need to be done:

      A B C D
    A   X X X
    B     X X
    C       X
    D

    That is, AB, AC, AD, BC, BD, and CD.

    Now we need to filter for transitive paths.  Let's start by considering
    only forward paths, i.e. paths where all matches are in the forward
    orientation.

    First, find all nearest-neighbor connections, paths A->B->C->D.
    This will be the set of possible paths, taken from the comparisons
    AB, BC, and CD (the outermost diagonal on the matrix above).

    We're left with the comparisons AC, AD, and BD.

    Now, for each of the possible A->B->C->D paths, verify that there
    is a direct connection from A->C, A->D, and B->D.  Voila!

    Things become more complicated when you have to take into account
    orientation; it involves counting orientations.  I'm not sure I
    understand it completely.
    """
    
    def __init__(self, w, t, *seqs):
        self.windowsize = w
        self.threshold = t
        self.seqs = [ ]
        self.cmps = {}

        for seq in seqs:
            self.add_sequence(seq)

    def add_sequence(self, seq):
        """
        """
        n = len(self.seqs)
        
        for i in range(0, n):
            cmp = paircomp.do_rolling_nxn_compare(self.seqs[i], seq,
                                                  self.windowsize,
                                                  self.threshold)
            d = self.cmps.get(i, {})
            d[n] = cmp
            self.cmps[i] = d

        self.seqs.append(seq)

    def make_paths(self, start_seq_i, pos):
        """
        Recursive function that builds a transitive path across
        A->B, B->C, C->D, ...
        """
        # exit?
        next_seq_i = start_seq_i + 1
        if next_seq_i == len(self.seqs):
            return [ [_PosO(pos, 1)] ]

        # retrieve matches for the given position
        cmp = self.cmps[start_seq_i][next_seq_i]
        matches = cmp[pos]

        # record each path element, along with its orientation, and recurse.
        paths = []
        for m in matches:
            if m.orientation == 1:
                bot = m.bot
                pos_i = _PosO(pos, 1)
            else:
                bot = m.bot - self.windowsize + 1
                pos_i = _PosO(pos, -1)
                
            l = self.make_paths(next_seq_i, bot)

            # insert this coordinate along with the appropriate orientation.
            for p in l:
                p.insert(0, pos_i)
                paths.append(p)

        return paths

    def filter(self):
        """
        Construct a list of all paths that are transitive across all
        of the analyses.
        """
        
        # build a list of all inner comparisons to check.
        pairs = []
        for j in range(2, len(self.seqs)):
            pairs += [ (i, i + j) for i in range(0, len(self.seqs) - j) ]

        # build paths from the outermost diagonal: AB, BC, CD, ...
        paths = []
        for i in range(0, len(self.seqs[0])):
            p = self.make_paths(0, i)
            paths += p

        # 'paths' now contains a list of our putative paths, those
        # that *could* be transitive, because they stretch from
        # A->B->C->D->...

        # Check to make sure that the path elements are shared
        # between all of the comparisons, e.g. A->C, A->D, C->D, ...
        # discard any paths that aren't completely shared.
        
        keep_paths = []
        for path in paths:
            keep_path = True
            for (from_seq, to_seq) in pairs:
                if not self.check_match(from_seq, to_seq, path):
                    keep_path = False
                    break
                
            # only keep those paths that match *all*.
            if keep_path:
                keep_paths.append(path)

#        keep_paths.sort()
        return keep_paths

    def check_match(self, from_seq, to_seq, path):
        """
        Check that the actual paircomp comparison between the sequences
        indexed by 'from_seq' and 'to_seq' contains the coordinates in
        the given path.
        """
        # retrieve the coordinates
        top_pos = path[from_seq]
        bot_pos = path[to_seq]

        # correct for orientation: two reverse matches == forward match.
        o = 1
        for i in range(from_seq, to_seq):
            o = o * path[i].orient

        # retrieve the actual comparison & matches
        cmp = self.cmps[from_seq][to_seq]
        matches = cmp[top_pos]

        # now check to see if this coordinate pair is contained in the
        # matches.
        for m in matches:
            bot = m.bot

            # correct for orientation
            if o == -1:
                bot = bot - self.windowsize + 1

            # if we find one, we're golden!
            if bot == bot_pos:
                return True

        return False

class _PosO:
    """
    Keep track of both a match position *and orientation*.
    """
    def __init__(self, pos, orient):
        self.pos = pos
        self.orient = orient

    def __repr__(self):
        sign = ' +-'
        sign = sign[self.orient]
            
        return "%s(%s)" % (self.pos, sign)

    def __int__(self):
        return self.pos

    def __cmp__(self, a):
        return cmp(self.pos, int(a))

def do_3way_regression(a, b, c, windowsize, threshold):
    """
    Run the current paircomp 3-way and return it in a similar
    path format for comparison purposes.
    """
    ab = paircomp.do_rolling_nxn_compare(a, b, windowsize, threshold)
    bc = paircomp.do_rolling_nxn_compare(b, c, windowsize, threshold)
    ac = paircomp.do_rolling_nxn_compare(a, c, windowsize, threshold)

    (ab, bc, ac) = paircomp.parser.filter_transitively(ab, bc, ac)

    paths = []
    for i in range(0, len(a)):
        matches = ab[i]
        for m in matches:
            bot = m.bot
            if m.orientation == -1:
                bot = bot - windowsize + 1
                apos = _PosO(i, -1)
            else:
                apos = _PosO(i, 1)

            matches2 = bc[bot]
            for m2 in matches2:
                bpos = _PosO(bot, m2.orientation)

                if m2.orientation == -1:
                    cpos = _PosO(m2.bot - windowsize + 1, 1)
                else:
                    cpos = _PosO(m2.bot, 1)

                paths += [ [apos, bpos, cpos] ]

    paths.sort()
    return paths

def test():
    run_all = 1

    seq = 'AAAAAAAAAA'
    seq2 = 'TTTTTTTTTT'

    if run_all:
        trials = (
            ( (seq, seq, seq), '[[0(+), 0(+), 0(+)]]' ),
            ( (seq, seq, seq2), '[[0(+), 0(-), 0(+)]]' ),
            ( (seq, seq2, seq), '[[0(-), 0(-), 0(+)]]' ),
            ( (seq, seq2, seq2), '[[0(-), 0(+), 0(+)]]' )
            )

        for (seqs, result) in trials:
            nway_cmp = NwayComparison(10, 1.0, *seqs)
            p = nway_cmp.filter()
            assert repr(p) == result

            p2 = do_3way_regression(seqs[0], seqs[1], seqs[2], 10, 1.0)
            assert p == p2

    #####

    seq = 'AAAAAAAAAA'
    seq2 = 'TTTTTTTTTT'
    
    if run_all:
        trials = (
            ( (seq, seq, seq, seq), '[[0(+), 0(+), 0(+), 0(+)]]' ),
            ( (seq, seq, seq, seq2), '[[0(+), 0(+), 0(-), 0(+)]]' ),
            ( (seq, seq, seq2, seq2), '[[0(+), 0(-), 0(+), 0(+)]]' ),
            ( (seq, seq2, seq2, seq2), '[[0(-), 0(+), 0(+), 0(+)]]' ),
            ( (seq, seq2, seq, seq), '[[0(-), 0(-), 0(+), 0(+)]]' ),
            ( (seq, seq, seq2, seq), '[[0(+), 0(-), 0(-), 0(+)]]' ),
            ( (seq, seq2, seq, seq2), '[[0(-), 0(-), 0(-), 0(+)]]' ),
            ( (seq, seq2, seq2, seq), '[[0(-), 0(+), 0(-), 0(+)]]' ),
            )

        for (seqs, result) in trials:
            nway_cmp = NwayComparison(10, 1.0, *seqs)
            assert repr(nway_cmp.filter()) == result

    #####

    if run_all:
        seq = 'AAAAAAAAAA'              # 10
        seq2 = 'TTTTTTTTTT'             # 10

        nway_cmp = NwayComparison(10, 1.0, seq, seq2, seq, seq)

        p = nway_cmp.filter()
        assert repr(p) == '[[0(-), 0(-), 0(+), 0(+)]]'

    if run_all:
        seq = 'AAAAAAAAAA'              # 10
        seq2 = 'TTTTTTTTTT'             # 10
        
        nway_cmp = NwayComparison(10, 1.0, seq, seq2, seq, seq, seq2)
        
        p = nway_cmp.filter()
        assert repr(p) == '[[0(-), 0(-), 0(+), 0(-), 0(+)]]'

    if run_all:
        seq = 'AAAAAAAAAAA'             # 11

        nway_cmp = NwayComparison(10, 1.0, seq, seq, seq)

        p = nway_cmp.filter()
        assert repr(p) == '[[0(+), 0(+), 0(+)], [0(+), 0(+), 1(+)], [0(+), 1(+), 1(+)], [0(+), 1(+), 0(+)], [1(+), 1(+), 1(+)], [1(+), 1(+), 0(+)], [1(+), 0(+), 0(+)], [1(+), 0(+), 1(+)]]'

        p2 = do_3way_regression(seq, seq, seq, 10, 1.0)
        
        p.sort()
        p2.sort()
        
        assert p == p2

    if run_all:
        seq = 'AAAAAAAAAAA'             # 11
        seq2 = 'TTTTTTTTTTT'            # 11

        nway_cmp = NwayComparison(10, 1.0, seq, seq2, seq)

        p = nway_cmp.filter()
        assert repr(p) == '[[0(-), 1(-), 0(+)], [0(-), 1(-), 1(+)], [0(-), 0(-), 1(+)], [0(-), 0(-), 0(+)], [1(-), 0(-), 1(+)], [1(-), 0(-), 0(+)], [1(-), 1(-), 0(+)], [1(-), 1(-), 1(+)]]'

        p2 = do_3way_regression(seq, seq2, seq, 10, 1.0)
        
        p.sort()
        p2.sort()
        
        assert p == p2

if __name__ == '__main__':
    test()
