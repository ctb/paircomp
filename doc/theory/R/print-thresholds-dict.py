#! /usr/bin/env python
from threshold_str import threshold_str

thresholds = threshold_str.split('\n')
thresholds = [ i.split() for i in thresholds ]
thresholds = [ (float(b), float(c), float(d)) for (a, b, c, d) in thresholds ]
thresholds = [ (int(a), int(b), int(c)) for (a, b, c) in thresholds ]

ws_dict = {}
for (seqsize, windowsize, cutoff) in thresholds:
    d = ws_dict.get(windowsize, {})
    d[seqsize] = cutoff
    ws_dict[windowsize] = d

k = ws_dict.keys()
k.sort()

d = ws_dict.values()[0]
k2 = d.keys()
k2.sort()

print ws_dict
