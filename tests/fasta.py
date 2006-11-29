"""
Simple FASTA parsing library.

functions:

* seq = load_single(f, strict=1, force=None) -- load a single sequence from f
* dict = load(f, strict=1, force=None)       -- load all the sequences in f
* is_protein(seq), is_DNA(seq)               -- test for valid protein/DNA
* force_protein(seq), force_DNA(seq)         -- replace invalid ch with X/N
* reverse_complement(dna_seq), rc(dna_seq)   -- make RC of DNA sequence
* reverse(dna_seq), complement(dna_seq)      -- reverse/complement DNA seq
"""

import string, re, sys, cStringIO
from array import array
import types

__version__ = 1.0
__printdebug__ = 0

legal_dna = "ACGTNX-"
legal_protein = "ABCDEFGHIKLMNPQRSTUVWXYZ-"

complementTranslation = string.maketrans('ACTG', 'TGAC')

def load_single(f, strict=1, force=None):
    d = load(f, strict, force)
    if len(d) != 1:
        raise FastaException("there can only be one sequence in file loaded with load_single! file %s" % (f,))
    return d.values()[0]

#
# load
#

def load(f, strict=1, force=None):
    """
    Loads sequences in FASTA format from the given file object.
    Returns dict.

    The 'strict' option forces the file to be in a stricter FASTA format,
    i.e. it must contain ONLY sequences, and they must be ONLY
    DNA or Protein sequences.

    The 'force' option is by default 'None'; if it is set to 'DNA' or
    'protein', all sequences are forced into that sequence type.  'force'
    implies 'strict'.
    """
    # check parameters:
    if not force is None:
        if not force in ('DNA', 'protein'):
            raise FastaParameterException("'force' must be 'DNA' or 'protein', not '%s'" % (force,))

        if not strict:                  # 'force' implies 'strict'.
            strict = 1
    
    # initialize empty dictionary, to be filled with sequences
    d= {}

    # read name

    try:
        success = 0
        try:
            l = f.readline()
            success = 1
        except AttributeError:              # is it maybe a filename?
            if type(f) == types.StringType and len(f) < 500:
                try:
                    f = open(f)
                    l = f.readline()
                    success = 1
                except IOError:
                    pass

            if not success:             # ok, see if it's interpretable as
                                        # a string.
                f = cStringIO.StringIO(f)
                l = f.readline()
                success = 1
    except:
        raise
    #pass

    if not success:
        raise FastaParameterException("don't know how to deal with first argument; must be either a file handle, a filename, or a string in FASTA format.")

    # check for Mac EOL characters.  Note that this assumes that the first
    # line we're getting is representative of the file, which will usually
    # be the case, because 'readline' seems to pay attention to platform-
    # specific EOLs.
    #
    # If we do find a Mac EOL (^M, character # 13) then reconstruct the
    # file into a cStringIO object.  Unless people are dealing with > 20 MB
    # of sequence, that should work fine ;).
    #
    macEOL = '%c' % (13,)
    if macEOL in l:
        l = l + f.read()                # get all of the file
        
        l = l.replace(macEOL, "\n")     # replace EOL characters with \n
        
        f = cStringIO.StringIO(l)       # replace file with StringIO

    #
    # Read through the file and pick out the good lines.  Ignore
    # non-FASTA gobbledygook at the top.
    #

    if not strict:
        while l and '>' not in l:
            l = f.readline()

        if l and '>' in l:              # reset to '>' in first place.
            i = l.find('>')
            l = l[i:]

    if not l:
        raise FastaFormatException("Error! no beginning '>'")

    if l[0] != '>':
        lines = l
        lines = lines + f.readline() + f.readline() + f.readline() + f.readline()
        lines = lines[0:500]            # cutoff at (arbitrary) 500...
        errorMsg = """

ERROR!  This doesn't look like a fasta file...

The first several lines are:
--------------------
%s
...
--------------------
</pre>

""" % (l,)
        raise FastaFormatException(errorMsg)

    nextName = l
    sequence = ""

    while nextName:
        name = string.rstrip(nextName[1:])
        while 1:
            l = f.readline()
            if not l:                   # end of file
                nextName = None
                break
            elif l[0] == '>':           # end of sequence
                nextName = l
                break

            sequence = sequence + string.rstrip(l)
        # end while 1

        sequence = string.upper(sequence)
        sequence = string.join(string.split(sequence), '')

        if len(name) and len(sequence):
            if strict:
                if force is None:
                    if (is_protein(sequence) or is_DNA(sequence)):
                        pass
                    else:
                        raise FastaContentException("Error, sequence '%s' contains illegal characteres.  Cannot determine if it's DNA or protein..." % (name,))
                else:                   # force is something
                    if force == 'DNA' and not is_DNA(sequence):
                        sequence = force_DNA(sequence)
                    elif force == 'protein' and not is_protein(sequence):
                        sequence = force_protein(sequence)
                        #raise "err"
                
            else:                       # only way we get here is if error...
                raise FastaContentException("Error, sequence contains illegal characters.")

            d[name] = sequence
            
            sequence = ""
            
    # end while name

    return d
# end load

#
# is_protein
#

def is_protein(seq):
    """
    Returns 1 if it contains only legal values for a protein sequence.

    c.f.  http://www.ncbi.nlm.nih.gov/BLAST/fasta.html
    """

    for ch in seq:
        if ch not in legal_protein:
            return 0

    return 1
# is_protein

#
# is_DNA
#

def is_DNA(seq):
    """
    Returns 1 if it contains only legal values for a DNA sequence.

    c.f.  http://www.ncbi.nlm.nih.gov/BLAST/fasta.html
    """
    for ch in seq:
        if ch not in legal_dna:
            return 0

    return 1
# is_DNA

is_dna = is_DNA

#
# force_DNA
#

def force_DNA(seq):
    """
    Removes all whitespace & then replaces all non-ACGT characters with
    N.
    """
    seq = string.join(string.split(seq.upper()))

    seq2 = []
    for ch in seq:
        if ch not in legal_dna:
            ch = 'N'

        seq2.append(ch)

    return string.join(seq2, '')
# force_DNA

force_dna = force_DNA

#
# force_protein
#

def force_protein(seq):
    """
    Removes all whitespace & then replaces all non-protein characters with
    X.
    """
    seq = string.join(string.split(seq))

    seq2 = []
    for ch in seq:
        if ch not in legal_protein:
            ch = 'X'

        seq2.append(ch)

    return string.join(seq2, '')
# force_protein

#
# reverse_complement
#

def reverse_complement(s):
    """
    Build reverse complement of 's'.
    """
    s = string.upper(s)
    assert is_DNA(s), "Your sequence must be DNA!"

    r = reverse(s)
    rc = complement(r)

    return rc

rc = reverse_complement                 # alias 'rc' to 'reverse_complement'

#
# complement
#

def complement(s):
    """
    Return complement of 's'.
    """
    c = string.translate(s, complementTranslation)
    return c

#
# reverse
#

def reverse(s):
    """
    Return reverse of 's'.
    """
    r = array('c', s)
    r.reverse()
    r = string.join(r, '')

    return r

class FastaException(Exception):
    pass

class FastaParameterException(FastaException):
    pass

class FastaFormatException(FastaException):
    pass

class FastaContentException(FastaException):
    pass

#############################################################################



#############################################################################

def _assert():                          # check stuff
    pass
