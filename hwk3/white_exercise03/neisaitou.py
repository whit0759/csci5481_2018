###
# Christohper White
# CSCI 5481 - Fall 2018
# November 1, 2018
###

import numpy as np
import pandas as pd

class Sequence(object):
    """
    Class that creates a Seqence object
    """

    def __init__(self, seqID, strSeq):
        self._seq = strSeq
        self._id = str(seqID)
        self.edge = None

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

    def __init__(self, name=None):
        self._edges = []
        self._name = name

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

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

    def __init__(self, head, tail, length=None):
        self.head = head
        self.tail = tail
        self.length = length

    @property
    def head(self):
        return self._head

    @head.setter
    def head(self, new_head):
        self._head = new_head

    @property
    def tail(self):
        return self._tail

    @tail.setter
    def tail(self, new_tail):
        self._tail = new_tail

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, new_len):
        self._length = new_len

    @property
    def name(self):
        return 'Edge_'+ self.head.name + '_' + self.tail.name

    
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
    
    keys = [seq.name for seq in seqs]

    return dict(zip(keys,seqs))

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

    return pd.DataFrame(columns=labels, data=dist_matrix, index=labels)

def calc_qmat(dist_matrix):
    """
    Calculate the Q-matrix from the distance matrix and return the Q-matrix
    """

    qmat = pd.DataFrame(data=np.zeros(dist_matrix.shape),
            index=dist_matrix.index, columns=dist_matrix.columns)

    (nrow, ncol) = qmat.shape
    assert nrow==ncol,'Matrix is not square'

    # The distance matrix is symmetric. Therefore the sum of all elements in
    # row i equals the sum of all elements in row j iff i==j
    delta_sum = dist_matrix.sum(axis=0)
    
    # Using the row/column sums that were already calculated will save computation
    for i in qmat.columns:
        for j in qmat.index:
            if i != j:
                qmat.loc[i,j] = (nrow-2)*dist_matrix.loc[i,j] - delta_sum.loc[i] - delta_sum.loc[j]

    return qmat
    

def build_tree(fastafile):

    nodeslist = read_fasta(fastafile)
    edgelist = {}

    dist_matrix = build_matrix(nodeslist.values())

    count = len(dist_matrix)
    nodeNum = 1

    while count>1:

        """ Calculate (joining) the Q-matrix """
        qmatrix = calc_qmat(dist_matrix)

        """ Branch length calculation """
        indx1,indx2 = np.unravel_index(qmatrix.values.argmin(), qmatrix.shape)
        lbl1 = qmatrix.columns[indx1]
        node1 = nodeslist[lbl1]
        lbl2 = qmatrix.columns[indx2]
        node2 = nodeslist[lbl2]

        deltaindx1 = (dist_matrix.iloc[indx1,indx2]/2 + 1/(2*(count-2)) *
                    (np.sum(dist_matrix.iloc[indx1,:]) - np.sum(dist_matrix.iloc[indx2,:])))

        deltaindx2 = dist_matrix.iloc[indx1,indx2] - deltaindx1
        
        assert deltaindx1>=0 and deltaindx2>=0, 'An edge length is less than zero'

        # Create and append edges and nodes
        newNode = 'Node'+nodeNum
        nodeslist[newNode] = Node(name=newNode)
        newEdge1 = Edge(lbl1, newNode, length=deltaindx1)
        edgelist[newEdge1.name] = newEdge1
        newEdge2 = Edge(lbl2, newNode, length=deltaindx2)
        edgelist[newEdge2.name] = newEdge2

        """ Update distance matrix """
        new_node_frame = 1/2 * (dist_matrix[lbl1]+
                                dist_matrix[lbl2]-
                                dist_matrix[lbl1][lbl2])

        new_node_frame.name = newNode
        new_node_frame[newNode] = 0
        new_dist_matrix = dist_matrix.append(new_node_frame)
        new_dist_matrix[newNode] = new_node_frame
        
        new_dist_matrix.drop([lbl1, lbl2], axis=0, inplace=True)
        new_dist_matrix.drop([lbl1, lbl2], axis=1, inplace=True)

        """
        My thought here is that I will drop the indx1, indx2 columns and rows. Then
        append the new node in row 0 and column 0. Then recalculate.
        """

        count -= 1
        nodeNum += 1

    return nodeslist, edgelist