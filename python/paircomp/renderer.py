# See README.txt for license and copyright information.

"""
Renderers for seqcomp sequence comparisons.

Exports two functions:
   * render_comparison(comparison, (top_start, top_end),
                                   (bot_start, bot_end),
                                   match_renderer_fn, passthru)
                                   
   * render_pixelized_comparison(comparison, (top_start, top_end),
                                             (bot_start, bot_end),
                                             top_pixels, bot_pixels,
                                             pixel_renderer_fn, passthru)

'comparison' is a Seqcomp/ShoudanComparison object, top_start, top_end,
bot_start, and bot_end are limits on what to render, and top_pixels/bot_pixels
are the number of pixels in which to render it.

'match_renderer_fn' is a Python function that takes 6 arguments:
    passthru (identical to object passed into render_comparison)
    top_pos of match (integer)
    bot_pos of match (integer)
    orientation of match (+1/-1)
    length of match (integer)
    number of matching bp (integer)

'pixel_renderer_fn' is a Python function that takes 4 arguments:
    passthru (identical to object passed into render_pixelized_comparison)
    pixel # of top match (integer)
    pixel # of bot match (integer)
    number of matching bp (integer)

It is the responsibility of these two functions to actually do graphics
plotting; they are called once for each match in the comparison that
is within the specified limits.

(This module is a Python wrapper around the _seqcomp_renderer.so C/C++
extension module.)
"""

import _paircomp_renderer

#
# render_comparison
#

def render_comparison(comparison, top_limits, bot_limits,
                      top_pixels, bot_pixels, fn, passthru):

    cObj = comparison._comparisonObj
    top_s, top_e = top_limits
    bot_s, bot_e = bot_limits

    _paircomp_renderer.render_comparison(cObj,
                                        top_s, top_e,
                                        bot_s, bot_e,
                                        fn, passthru)

#
# render_pixelized_comparison
#

def render_pixelized_comparison(comparison, top_limits, bot_limits,
                                top_pixels, bot_pixels, fn, passthru):

    cObj = comparison._comparisonObj
    top_s, top_e = top_limits
    bot_s, bot_e = bot_limits

    _paircomp_renderer.render_pixelized_comparison(cObj,
                                                  top_s, top_e, top_pixels,
                                                  bot_s, bot_e, bot_pixels,
                                                  fn, passthru)
