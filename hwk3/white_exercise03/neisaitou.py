###
# Christohper White
# CSCI 5481 - Fall 2018
# November 1, 2018
###

import numpy as np

class Sequence(object):
    """
    Class that creates a Seqence object
    """

    def __init__(self, seqID, strSeq):
        self._seq = strSeq
        self._id = seqID
        self._edge = None

    @property
    def name(self):
        return str(self._id)
    
    @property
    def seq(self):
        return self._seq

    @property
    def edge(self):
        return self._edge

    @edge.setter
    def edge(self, new_edge):
        self._edge = new_edge

    def __str__(self):
        return "Sequence {0.name!s} begins {1!s}...".format(self, self.seq[:10])

class Node(object):
    """
    Class for the internal nodes of the tree
    """

    def __init__(self):
        self._edges = []

    @property
    def edges(self):
        return self._edges

    def append(self, edge):
        self._edges.append(edge)
        edge.append(self)

    def remove(self, edge):
        edge.remove(self)
        try:
            self._edges.remove(self)
        except ValueError:
            print('The edge {} is not in node {}'.format(edge, self))

class Edge(object):
    """
    Class for the edges of the tree
    """

    def __init__(self):
        self._nodes = []
        self._length = None

    @property
    def nodes(self):
        return self._nodes

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, new_len):
        self._length = new_len

    def append(self, node):
        assert len(self._nodes) < 2, 'Cannot add any more nodes to {}.'.format(self)
        
        self._nodes.append(node)
        
    def remove(self, node):

        try:
            self._nodes.remove(node)
        except ValueError:
            print('The node {} is not in edge {}.'.format(node, self))

    

def read_fasta(file):
    """
    Function to read in a FASTA file and return a list of
    sequences in the file.
    """

    seqs = []

    with open(file, 'r') as f:
        line = f.readline()
        while line:
            if line[0]=='>':
                # Line is a sequence identifier
                seqid = line[1:].strip()
                
                # Read next line as the sequence
                line = f.readline()
                seq = line.strip()

                seqs.append(Sequence(seqid,seq))

                line = f.readline()
    
    return seqs

def calc_dist(seq1, seq2):
    """
    Calculate the similarity of two sequences and return the distance.
    """

    assert(len(seq1)==len(seq2))
    l = len(seq1)

    count = 0
    for n in range(l):
        if seq1[n].casefold()==seq2[n].casefold():
            # If the character at position 'n' is the same in both sequences,
            # increment the counter
            count += 1

    dist = l - count

    return dist


def build_matrix(seqlist):
    """
    Build the initial distance matrix for the sequences in seqlist. Return the matrix.
    """

    l = len(seqlist)

    dist_matrix = np.zeros((l,l))
    labels = [ x.name for x in seqlist ]

    for n in range(l-1):
        for m in range(n+1,l):
            dist_matrix[n,m]=dist_matrix[m,n]=calc_dist(seqlist[n].seq,seqlist[m].seq)

    return labels, dist_matrix

def calc_qmat(dist_matrix):
    """
    Calculate the Q-matrix from the distance matrix and return the Q-matrix
    """

    qmat = np.zeros_like(dist_matrix)

    (nrow, ncol) = np.shape(dist_matrix)
    assert(nrow==ncol)

    # The distance matrix is symmetric. Therefore the sum of all elements in
    # row i equals the sum of all elements in row j iff i==j
    delta_sum = np.sum(dist_matrix, axis=0)

    # Using the row/column sums that were already calculated will save computation
    for i in range(nrow):
        for j in range(ncol):
            qmat[i,j] = (nrow-2)*dist_matrix[i,j] - delta_sum[i] - delta_sum[j]

    return qmat
    

def build_tree(fastafile):

    seqlist = read_fasta(fastafile)
    nodeslist = []
    edgelist = []

    labels, dist_matrix = build_matrix(seqlist)

    while len(labels)>1:

        """ Calculate (joining) the Q-matrix """
        qmatrix = calc_qmat(dist_matrix)

        """ Branch length calculation """
        indx1,indx2 = np.min(qmatrix)
        
        deltaindx1 = (dist_matrix[indx1,indx2]/2 + 1/(2*(len(labels)-2)) *
                    (np.sum(dist_matrix[indx1,:]) - np.sum(dist_matrix[indx2,:])))

        deltaindx2 = dist_matrix[indx1,indx2] - deltaindx1
        
        if deltaindx1 < 0 or deltaindx2 < 0:
            deltaindx1 = (dist_matrix[indx2,indx1]/2 + 1/(2*(len(labels)-2)) *
                        (np.sum(dist_matrix[indx2,:]) - np.sum(dist_matrix[indx1,:])))

            deltaindx2 = dist_matrix[indx2,indx1] - deltaindx1
    
        assert deltaindx1>=0 and deltaindx2>=0, 'An edge length is less than zero'

        ## Create and append edges and nodes

        """ Update distance matrix """
        old_dist_matrix = dist_matrix

        """
        My thought here is that I will drop the indx1, indx2 columns and rows. Then
        append the new node in row 0 and column 0. Then recalculate.
        """

