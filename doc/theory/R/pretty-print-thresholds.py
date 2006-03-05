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

print '''
<h3>Cutoffs for random background</h3>
The thresholds/windowsizes below are the thresholds at which you
have a 5% chance of seeing a match, in random DNA of the length given,
to a window of the given windowsize.
<p>
In the table below, windowsizes are on the left, and
sequence sizes are along the top.
<p>
<table border=1>
'''

print '<tr><td></td>',
for seqsize in k2:
    print '<td><b>%d</b></td>' % (seqsize,),
print '</tr>'

for ws in k:
    print '<tr><td><b>%d</b></td>' % (ws,),
    d = ws_dict[ws]
    for seqsize in k2:
        cutoff = d[seqsize]
        percent = float(cutoff) / float(ws) * 100
        print '<td align=center>%d%%</td>' % (percent),

    print '</tr>'

print '</table>'

#print 'size\tws\tcutoff\tpercent'
#for seqsize, windowsize, cutoff in thresholds:
#    print "%-6d\t%-2d\t%-2d\t%-3d%%" % (seqsize, windowsize, cutoff,
#                                     float(cutoff) / float(windowsize) * 100)
