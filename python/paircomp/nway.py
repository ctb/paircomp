import _paircomp_parser

class NwayComparison:
    def __init__(self, windowsize, threshold, *seqs):
        self.windowsize = windowsize
        self.threshold = threshold

        self._nway = _paircomp_parser.create_nway(windowsize, threshold)
        for seq in seqs:
            self.add_sequence(seq)

    def add_sequence(self, seq):
        _paircomp_parser.add_sequence_to_nway(self._nway, seq)

    def get_filtered_paths(self):
        paths = _paircomp_parser.get_nway_filtered_paths(self._nway)

        l = []
        for path in paths:
            j = [ apply(_PositionAndOrientation, element) for element in path ]
            l.append(j)
            
        return l

class _PositionAndOrientation:
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
