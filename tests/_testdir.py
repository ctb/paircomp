"""
Inform the test scripts of the path to paircomp's Python ext module.

This is cross platform.
"""

import sys
import os.path
import distutils.util
import platform

#
# get the current directory from __file__
#

testdir = os.path.abspath(os.path.dirname(__file__))

#
# now build paircomp's python build/lib directory path
#

platform_str = distutils.util.get_platform()
version = ".".join(platform.python_version_tuple()[:2])
paircomp_build_dir = testdir + '/../python/build/lib.%s-%s' % (platform_str,
                                                             version)
paircomp_build_dir = os.path.abspath(paircomp_build_dir)

#
# put it in the path
#

sys.path.insert(0, paircomp_build_dir)

import paircomp
assert paircomp.__version__ == "1.0"

def paircomp_version_message():
    print """
Testing version %s of the paircomp extension module, loaded from

\t%s
""" % (paircomp.__version__, paircomp_build_dir,)


