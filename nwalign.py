# Generate the alignment of two sequences using the Needleman-Wunsch algorithm
#
# Author: Christopher White
# Date: September 29, 2018
#


import sys, os
import argparse
import numpy as np

GAP = -2
MATCH = 1
MISMATCH = -3

def make_arg_parser():
    parser = argparse.ArgumentParser(prog='nwalign.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--match",
                        default=argparse.SUPPRESS,
                        required=False,
                        help="Matches file for anchored alignment")
    parser.add_argument("-o","--output",
                      default=None,
                      required=False,
                      help="Path to output file")
    parser.add_argument("-V","--verbose",
                      action="store_true",
                      help="Verbose output")
    return parser

def gen_matrix(size=()):
    """
    Return the initial seeded alignment matrix
    """
    nx, ny = size

    align_mat = np.zeros(size)
    horz, vert = np.meshgrid(np.arange(ny), np.arange(nx), sparse=True)
    
    align_mat[:1,:] = horz
    align_mat[:,:1] = vert
    
    align_mat = align_mat*-2
    
    return align_mat

def read_seq(seqfile):
    """
    Returns the sequence in seqfile. The format is expected to be a header line
    by a sequence wrapped onto multiple lines.
    """

    with open(seqfile) as f:
        f.readline()
        seq = f.read()

    seq = seq.replace('\n', '')

    return seq

def read_match(matchfile):
    """
    Returns the sequence anchors read from a match file.
    """
    
    anchor_mat = np.loadtxt(matchfile)
    
    human_match = anchor_mat[:,:2]
    fly_match = anchor_mat[:,2:]
    
    return fly_match, human_match

def align_matrix(seq1, seq2):
    """
    Returns the alignment matrix using Needleman-Wunsch
    """

    nh = len(seq1)
    nv = len(seq2)

    nwmat = gen_matrix(size=(nv+1,nh+1))

    gap = GAP

    for i in range(nv):
        """ Start with operating on each row """
        for j in range(nh):
            """ Step through each entry of each row """
    
            if seq2[i]==seq1[j]:
                match = MATCH
            else:
                match = MISMATCH
                
            nwmat[i+1, j+1] = max(nwmat[i,j]+match, nwmat[i+1,j]+gap, 
                 nwmat[i,j+1]+gap)

    return nwmat

def find_alignments(nwmat,seq1,seq2):
    """
    Returns the alignment score and the optimal alignments of the sequences
    """
    
    # The score is the entry of the last row and last column
    score = nwmat[-1,-1]
    
 
    newseq1 = ''
    newseq2 = ''
    
    i = len(seq2)
    j = len(seq1)
     
    gap = GAP
    
    while i>0 or j>0:
        # Step through the columns and rows in reverse order
        # Recalculate the route backwards
        if seq1[j-1] == seq2[i-1]:
            match = MATCH
        else:
            match = MISMATCH
        
        # While not in the first row or first column, compare the three 
        # adjacent values to the current value. If in the first row or first
        # column, then prepend '_' (gaps) until we reach (0,0)
        if i>0 and j>0:
            if nwmat[i-1, j-1] == (nwmat[i,j]-match):
                i-=1
                j-=1
                newseq2 = seq2[i] + newseq2
                newseq1 = seq1[j] + newseq1
            elif nwmat[i-1,j] == (nwmat[i,j]-gap):
                i-=1
                newseq2 = '_' + newseq2
            elif nwmat[i,j-1] == (nwmat[i,j]-gap):
                j-=1
                newseq1 = '_' + newseq1
        elif i>0:
            i-=1
            newseq2 = '_' + newseq2
        elif j>0:
            j-=1
            newseq1 = '_' + newseq1
    
    return score, newseq1, newseq2
