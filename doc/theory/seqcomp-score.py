"""
A simple scoring system for seqcomp-style pairwise analyses.

CTB 7/2003.
"""

import math

def log2(n):
    "log base 2."
    return math.log(n) / math.log(2)

def log_match(l):
    assert l >= 0 and l < 1, "invalid value for lambda"
    P_match = (3. * l + 1.) / 16.
    return log2(P_match)

def log_mismatch(l):
    assert l >= 0 and l < 1, "invalid value for lambda"
    P_mismatch = (1. - l) / 16.
    return log2(P_mismatch)

def score(M, N, l=.95):
    """
    Return a score for M matches in N bases with parameter 0 <= l < 1.
    l defaults to 0.95.
    """
    
    M = int(M)
    N = int(N)
    l = float(l)
    
    sum_random = -4 * N
    sum_pair = M * log_match(l) + (N - M) * log_mismatch(l)

    return sum_pair - sum_random

#
# Run from the command line.
#

if __name__ == '__main__':
    print """
To see some scores, try the following:

---
>> print score(15, 20)
>> print score(15, 20, 0)
---

where '15' is the number of matches out of a window of '20', and the third
(optional) parameter is the lambda.  By default, lambda is set to 0.95.
"""
