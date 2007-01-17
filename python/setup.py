from distutils.core import setup, Extension

# the c++ parser extension module (needs to be linked in with libpaircomp.a...)
parser_mod = Extension("paircomp._paircomp_parser",
                       ["c++-ext/_paircomp_parser.cc",],
                       include_dirs=['../lib',],
                       library_dirs=['../lib',],
                       libraries=['paircomplib', 'stdc++'],
                       depends=['../lib/libpaircomplib.a',])

# the c++ algorithm extension module (needs to be linked in with libpaircomp.a...)
algorithms_mod = Extension("paircomp._paircomp_algorithms",
                           ["c++-ext/_paircomp_algorithms.cc",],
                           include_dirs=['../lib',],
                           library_dirs=['../lib',],
                           libraries=['paircomplib', 'stdc++'],
                           depends=['../lib/libpaircomplib.a',])

# the c++ renderer extension module
#renderer_mod = Extension("paircomp._paircomp_renderer",
#                         ["c++-ext/_paircomp_renderer.cc"],
#                         include_dirs=['../lib',],
#                         library_dirs=['../lib',],
#                         libraries=['paircomplib', 'stdc++'],
#                         depends=['../lib/libpaircomplib.a',])

# python modules
package = 'paircomp'

setup(name = "paircomp", version = "1.0",
      description = 'paircomp sequence comparison library',
      author = 'C. Titus Brown',
      author_email = 'titus@caltech.edu',
      url = 'http://family.caltech.edu/',
      packages = [package,],
      ext_modules = [parser_mod, algorithms_mod])
