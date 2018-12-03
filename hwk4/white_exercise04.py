"""
CSCI 5481
Homework 4

Christopher White
November 20, 2018
"""

import sys
import numpy as np 
import pandas as pd 

class Sequence(object):
    """
    Class that creates a Seqence object
    """

    def __init__(self, seqID, strSeq):
        self._seq = strSeq
        self._id = seqID

    @property
    def seq(self):
        return self._seq

    @property
    def name(self):
        return self._id

    @property
    def data_frame(self):
        return pd.DataFrame(data=list(self.seq), columns=[self.name], dtype=np.unicode_)

    def __str__(self):
        return "Sequence {0.name!s} begins {1!s}...".format(self, self.seq[:10])

def read_fasta(file):
    """
    Function to read in a FASTA file and return a list of
    sequences in the file.
    """

    seqs = []
    count = 0

    with open(file, 'r') as f:
        line = f.readline()
        while line:
            if line[0]=='>':
                # Line is a sequence identifier
                seqid = line[1:].strip()
                
                # Read next line as the sequence
                line = f.readline()
                seq = line.strip()

                count += 1
                seqs.append(Sequence(seqid,seq))

                line = f.readline()
    
    keys = [seq.name for seq in seqs]

    return dict(zip(keys,seqs))

def print_regions(seqmin, seqmax, ind):
    """
    Prints the variable regions to a file given a min value, max value
    and an array of positions.
    The positions array are all positions that are in the variable
    regions.
    """

    # Find the region boundaries by iterating a cursor across the array.
    # When the cursor is greater than the last position + 1, then the last
    # position is the end of a region and the cursor is the start of the
    # next region.
    last = ind.min()
    starts = [last]
    ends = []
    for ind in ind:
        if last+1<ind:
            ends.append(last)
            starts.append(ind)
        
        last = ind

    ends.append(last)


    # Write the regions to file. If the region begins at the beginning of
    # the sequence or ends at the end of the sequence, then don't write it.
    # If the region is smaller than 3 positions, then don't write it.
    with open('hwk4_prob3_regions.csv', 'w') as f:
        for start,end in zip(starts,ends):
            if start==seqmin or end==seqmax or (end-start < 3):
                continue

            f.write('{0}, {1}\n'.format(start, end))


if __name__ == '__main__':

    import numpy as np
    import matplotlib as mpl 
    import matplotlib.pyplot as plt


    if len(sys.argv)==2:
        fna_file = sys.argv[1]
    else:
        fna_file = 'Homework_4_seqs.fna'

    # Read in sequences from the file
    seqs = read_fasta(fna_file)

    num_seqs = len(seqs)

    # Create a dataframe that has the sequences as keys
    seqdf = pd.concat([seqs[key].data_frame for key in seqs.keys()],
                        axis='columns')

    print(seqdf.head().iloc[:,:5])

    # Count the number of each nucleic acid in each position.
    # Create a new dataframe that is keyed by position with the number
    # of acids of each type
    acidcount = pd.DataFrame({'A': seqdf.where(seqdf=='A').T.count()})
    acidcount['C'] = seqdf.where(seqdf=='C').T.count()
    acidcount['G'] = seqdf.where(seqdf=='G').T.count()
    acidcount['T'] = seqdf.where(seqdf=='T').T.count()
    acidcount['-'] = seqdf.where(seqdf=='-').T.count()

    # Find the most common type at each position in the sequence
    freq = acidcount[['A','C','G','T']].max(axis=1).div(num_seqs)

    # Smooth the data
    data = freq.rolling(32, center=True, min_periods=8).mean()

    data.to_csv('./hwk4_prob1_identity_values.csv')

    # Cube the data to make it more uniformly distributed
    #   The histogram of the data looked skewed toward 1 so
    #   I decided to find the threshold of a more uniformly
    #   distributed version of the data by taking data cubed.
    datacube = data.pow(3)

    # Find the mean and std of the cube of the data
    dcube_mean = datacube.mean()
    dcube_std = datacube.std()

    # Calculate the threshold for determining variable regions.
    #   Since I found the threshold of the cube of the data, I
    #   need to take the cube root for finding the variable regions
    #   from the original data.
    data_threshold = dcube_mean ** (1/3)

    # Plot the variation versus position
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_title('Variability versus Position for Given Sequences')
    ax.set_ylabel('Conserved Fraction')
    ax.set_xlabel('Sequence Position')

    data.plot(kind='line', ax=ax)
    plt.hlines(data_threshold, 0, len(data)-1)

    plt.draw()

    fig.savefig('hwk4_prob2_graph.pdf')


    # Print the variable regions to a CSV file
    print_regions(data.index.min(), data.index.max(), 
                    data.loc[data<data_threshold].index)

