# Generate the alignment of two sequences using the Needleman-Wunsch algorithm
#
# Author: Christopher White
# Date: September 29, 2018
#


import sys
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
                        default=None,
                        required=False,
                        help="Match file for anchored alignment")
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
    
    # INITIAL ALIGNMENT MATRIX
    # _ a b c ...
    # a 0 -2 -4 ...
    # b -2 0 0 ...
    # c -4 0 0 ...
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
    
    anchor_mat = np.loadtxt(matchfile, dtype='i4')
    
    # First two columns are for human sequence
    human_match = anchor_mat[:,:2]
    
    # Second two columns are for fly sequence
    fly_match = anchor_mat[:,2:]
    
    # Assert that the anchor pairs are equal in dimension
    assert(np.shape(human_match)==np.shape(fly_match))
    
    # Subtract 1 to account for 0 based indexing
    return human_match-1, fly_match-1

def align_matrix(seq1, seq2):
    """
    Returns the alignment matrix using Needleman-Wunsch
    """

    nh = len(seq1)
    nv = len(seq2)

    # Create initial matrix
    nwmat = gen_matrix(size=(nv+1,nh+1))

    gap = GAP

    for i in range(nv):
        """ Start with operating on each row """
        for j in range(nh):
            """ Step through each entry of each row """
    
            if seq2[i]==seq1[j]:
                # If the entries match, set match to MATCH
                match = MATCH
            else:
                match = MISMATCH
            
            # The matrix cell is the max of three conditions, either a match
            # or mistmatch plus gap penalties.
            nwmat[i+1, j+1] = max(nwmat[i,j]+match, nwmat[i+1,j]+gap, 
                 nwmat[i,j+1]+gap)

    return nwmat

def find_alignments(seq1,seq2):
    """
    Returns the alignment score and the optimal alignments of the sequences
    """
    nwmat = align_matrix(seq1, seq2)
    
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
            # If a match condition exists decrement i and j
            if nwmat[i-1, j-1] == (nwmat[i,j]-match):
                i-=1
                j-=1
                newseq2 = seq2[i] + newseq2
                newseq1 = seq1[j] + newseq1
            
            # If the current cell resulted from the cell above, insert a gap
            # in seq2
            elif nwmat[i-1,j] == (nwmat[i,j]-gap):
                i-=1
                newseq1 = '_' + newseq1
                newseq2 = seq2[i] + newseq2
            
            # If The current cell resulted from the cell to the left, insert a
            # gap in seq1
            elif nwmat[i,j-1] == (nwmat[i,j]-gap):
                j-=1
                newseq1 = seq1[j] + newseq1
                newseq2 = '_' + newseq2
        elif i>0:
            i-=1
            newseq1 = '_' + newseq1
            newseq2 = seq2[i] + newseq2
        elif j>0:
            j-=1
            newseq1 = seq1[j] + newseq1
            newseq2 = '_' + newseq2
    
    return (score, newseq1, newseq2)

def anchor_align(human_seq, fly_seq, human_match, fly_match):
    """
    Returns aligned sequences and anchored sequences for a given set of human
    and fly sequences with anchored regions given by 'match' arrays.
    """
    
    # The matches are given in start-stop pairs. The first aligment is done 
    # from the beginning to the first start. Then, alignments are between each
    # stop and the next start. I am assuming that an equal number of pairs are
    # given.
    
    # Assert that the max or min numbers in the match arrays don't exceed the
    # sequence bounds
    assert(np.max(human_match)<=len(human_seq))
    assert(np.max(fly_match)<=len(fly_seq))
    assert(np.min(human_match)>=0)
    assert(np.min(fly_match)>=0)
       
    # For all the sections between a 'stop' and a 'start' align the sections
    score = 0.0
    new_human = human_seq[:int(human_match[0,0])]
    new_fly = fly_seq[:int(fly_match[0,0])]

    num_matches = np.shape(human_match)[0]

    for n in range(num_matches):

        # Create the indicices
        hstart = int(human_match[n,0])
        hstop = int(human_match[n,1])
        fstart = int(fly_match[n,0])
        fstop = int(fly_match[n,1])
        
        # Align the matched sequence
        (tmpscore, _, _) = find_alignments(human_seq[hstart:hstop],
                                                     fly_seq[fstart:fstop])

        # Append the matched results
        score += tmpscore
        new_human += human_seq[hstart:hstop]
        new_fly += fly_seq[fstart:fstop]
        
        # Align the sequence between stop and next start
        if n<num_matches-1:
            hnext = int(human_match[n+1,0])
            fnext = int(fly_match[n+1,0])
            (tmpscore, tmphseq, tmpfseq) = find_alignments(human_seq[hstop:hnext],
                                                         fly_seq[fstop:fnext])
        
            # Append the matched results
            score += tmpscore
            new_human += tmphseq
            new_fly += tmpfseq
        
        elif n==num_matches-1:
            # Append the final sequences
            new_human += human_seq[hstop:]
            new_fly += fly_seq[fstop:]        
    
    return (score, new_human, new_fly)

def print_alignments(file, results):
    """
    Print the alignment scores and sequences
    """
    import textwrap
    
    score, seq1, seq2 = results
    
    print('Score: {}'.format(score), file=file)
        
    seq1wrap = textwrap.wrap(seq1, width=50)
    seq2wrap = textwrap.wrap(seq2, width=50)
    
    print('\nHUMAN SEQUENCE', file=file)
    print('*'*20, file=file)
    
    for line1 in seq1wrap:
        print('{:<}'.format(line1), file=file)
        
    print('\nFLY SEQUENCE', file=file)
    print('*'*20, file=file)
    
    for line2 in seq2wrap:    
        print('{:<}'.format(line2), file=file)
        
        

def main(query, ref, match=None, output=None):
    """
    Run the alignments for the given input files and if match or output files
    are given, process those as well.
    """
    
    # Read the sequences
    human_seq = read_seq(query)
    fly_seq = read_seq(ref)
    
    if match:
        # If a match file is given use it to do anchored alignment
        human_match, fly_match = read_match(match)
        results = anchor_align(human_seq, fly_seq, human_match,
                                          fly_match)
    else:
        # If no match file, align the whole sequences
        results = find_alignments(human_seq, fly_seq)
        
    if output:
        # If an output file is given, write to it
        with open(output,mode='w') as f:
            print_alignments(f, results)
            
    else:
        # If no output file, print to stdout
        print_alignments(sys.stdout, results)
        
    

if __name__=="__main__":
    
    parser = make_arg_parser()
    args = parser.parse_args()
    
    main(args.query, args.ref, match=args.match, output=args.output)