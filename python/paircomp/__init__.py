# Base package imports for the 'paircomp' Python package.
#
# See README.txt for license and copyright information.

__version__ = "1.0"

from parser import parse_paircomp, parse_seqcomp, filter_transitively, \
     build_transitive
#from renderer import render_comparison, render_pixelized_comparison
from algorithms import do_simple_nxn_compare, do_rolling_nxn_compare, \
     do_hashed_n_compare
