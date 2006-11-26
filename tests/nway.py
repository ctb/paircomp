import _testdir
import paircomp

class NwayComparison:
    def __init__(self, w, t, seq1):
        self.windowsize = w
        self.threshold = t
        self.seqs = [ seq1 ]
        self.cmps = {}

    def add_sequence(self, seq):
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
        if (start_seq_i + 1) == len(self.seqs):
            return [ [pos] ]

        cmp = self.cmps[start_seq_i][start_seq_i + 1]
        matches = cmp[pos]

        paths = []
        for m in matches:
            if m.orientation == 1:      # @CTB forward only
                l = self.make_paths(start_seq_i + 1, m.bot)
                for p in l:
                    p.insert(0, pos)
                    paths.append(p)

        return paths

    def filter(self):
        # build a list of all inner comparisons to check.
        pairs = []
        for j in range(2, len(self.seqs)):
            pairs += [ (i, i + j) for i in range(0, len(self.seqs) - j) ]

        # build paths from the outermost diagonal: AB, BC, CD, ...
        paths = []
        for i in range(0, len(self.seqs[0])):
            p = self.make_paths(0, i)
            paths += p

        # for each path, check each of the inner comparisons, too.
        keep_paths = []
        for path in paths:
            keep_path = True
            for (from_seq, to_seq) in pairs:
                top_pos = path[from_seq]
                bot_pos = path[to_seq]
                if not self.check_match(from_seq, to_seq, top_pos, bot_pos):
                    keep_path = False
                    break

            if keep_path:
                keep_paths.append(path)

        return keep_paths

    def check_match(self, from_seq, to_seq, top_pos, bot_pos):
        cmp = self.cmps[from_seq][to_seq]
        matches = cmp[top_pos]
        for m in matches:
            if m.orientation == 1 and m.bot == bot_pos:
                return True

class NwayPathPair:
    def __init__(self, top, bot, orient):
        self.top = top
        self.bot = bot
        self.orient = orient
        
class NwayFullPath:
    def __init__(self):
        self.d = {}

    def add_path(self, from_seq, to_seq, path_pair):
        e = self.d.get(from_seq, {})
        e[to_seq] = path_pair
        self.d[from_seq] = e

def test():
    seq = 'AAAAAAAAAAA'
    
    nway_cmp = NwayComparison(10, 1.0, seq)
    nway_cmp.add_sequence(seq)
    nway_cmp.add_sequence(seq)

    p = nway_cmp.filter()
    assert p == [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 1], [1, 1, 0], [1, 0, 0], [1, 0, 1]]
    
    assert nway_cmp.check_match(0, 1, 0, 0)
    assert nway_cmp.check_match(0, 1, 1, 1)
    assert nway_cmp.check_match(0, 1, 0, 1)
    assert nway_cmp.check_match(0, 1, 1, 0)

    nway_cmp.add_sequence(seq)
    nway_cmp.filter()

#    seq2 = 'TTTTTTTTTT'
    
#    nway_cmp = NwayComparison(10, 1.0, seq)
#    nway_cmp.add_sequence(seq2)
#    nway_cmp.add_sequence(seq)

#    p = nway_cmp.filter()
#    print p

if __name__ == '__main__':
    test()
