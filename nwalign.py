# Generate the alignment of two sequences using the Needleman-Wunsch algorithm
#
# Author: Christopher White
# Date: September 29, 2018
#


import sys, os
import parser
import numpy as np

def make_arg_parser():
    parser = argparse.ArgumentParser(prog='exercise02.py',
                          # version="%prog 1.0",
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-o","--output",
                      default=None,
                      required=False,
                      help="Path to output file [required]")
    parser.add_argument("-V","--verbose",
                      action="store_true",
                      help="Verbose output")
    return parser

def gen_matrix(size=()):
    """
    Return the initial alignment matrix
    """
    nx, ny = size

    align_mat = np.zeros(size, dtype=int16)
    align_mat[0,:], align_mat[:,0] = np.meshgrid(0:nx, 0:ny, sparse=True)
    align_mat = align_mat*-2
    return align_mat

def read_seq(seqfile):
    """
    Returns the sequence in seqfile
    """

    with open(seqfile) as f:
        pass

    return seq

def align_matrix(seq1, seq2):
    """
    Returns the alignment matrix using Needleman-Wunsch
    """

    nx = length(seq1)
    ny = length(seq2)

    nwmat = gen_matrix(size=(nx,ny))

    for i in range(nx):
        for j in range(ny):
            pass

    return nwmat

def find_alignments(nwmat):
    pass
    
