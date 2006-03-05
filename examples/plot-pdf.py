#! /usr/bin/env python
"""
Generates a PDF dot-plot from a given comparison, using the
seqcomp_utils module parser & renderer functions.
"""

comparison_file = '../rawdata/apo_bec/sp-lv.txt'
comparison_type = 'shoudan'             # shoudan or seqcomp
windowsize = 10
top_len = 44522                         # must be specified for Shoudan-style

# output file
pdf_file = 'test.pdf'

# plot from where to where?
top_limits = (0, 10000)
bot_limits = (0, 10000)

# at what resolution? (total pixels across which to plot)
x_pixels = 500
y_pixels = 500

###########################################################################3
import sys
sys.path.append('../')
import seqcomp_utils
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

width, height = letter

x_start = 1.5*inch
x_height = width - 3*inch

y_start = 1.5*inch
y_height = height - 3*inch

#
# Read stuff in.
#

r = open(comparison_file).read()
if comparison_type == 'shoudan':
    comparison = seqcomp_utils.ShoudanComparison(r, windowsize, top_len)

else:
    comparison = seqcomp_utils.SeqcompComparison(r)

# Initialize the canvas.

c = canvas.Canvas(pdf_file, pagesize=letter)

# draw a box around the entire canvas.
    
c.line(x_start, y_start, x_start, y_height + y_start)
c.line(x_start, y_height + y_start, x_start + x_height, y_height + y_start)
c.line(x_start + x_height, y_height + y_start, x_start + x_height, y_start)
c.line(x_start + x_height, y_start, x_start, y_start)

# calculate the scaling parameters

x_scale = float(x_height) / float(x_pixels)
y_scale = float(y_height) / float(y_pixels)

# create the function that will actually do the plotting:

def f(canvas, x, y, m):
    x_1 = x * x_scale + x_start
    x_2 = x_1 + x_scale
    
    y_1 = y * y_scale + y_start
    y_2 = y_1 + y_scale

    canvas.line(x_1, y_1, x_2, y_2)

# plot

print 'rendering'
seqcomp_utils.render_pixelized_comparison(comparison,
                                          top_limits, bot_limits,
                                          x_pixels, y_pixels,
                                          f, c)

# print to page & save

c.showPage()
c.save()
