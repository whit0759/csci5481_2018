###
# Christohper White
# CSCI 5481 - Fall 2018
# November 1, 2018
###

import sys
import numpy as np
import pandas as pd

class Node(object):
    """
    Class for the internal nodes of the tree
    """

    def __init__(self, name=None, index=None):
        self._edges = []
        self._name = name
        self._index = index
        self._children = []
        self._parent = None

    @property
    def children(self):
        return self._children
    
    @property
    def parent(self):
        return self._parent

    def add_child(self, new_child):
        self._children.append(new_child)

    def remove_child(self, child):
        self._children.remove(child)
    
    @parent.setter
    def parent(self, new_parent):
        self._parent = new_parent

    @property
    def index(self):
        return self._index

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def edges(self):
        return self._edges

    def append_edge(self, edge):
        self._edges.append(edge)
        
    def remove_edge(self, edge):
        try:
            self._edges.remove(self)
        except ValueError:
            print('The edge {} is not in node {}'.format(edge, self))

class Sequence(Node):
    """
    Class that creates a Seqence object
    """

    def __init__(self, seqID, strSeq, index=None):
        super().__init__(name=seqID, index=index)
        self._seq = strSeq
        self._children = None

    @property
    def seq(self):
        return self._seq

    @property
    def children(self):
        return self._children

    def add_child(self, new_child):
        if new_child is not None:
            raise AttributeError('Cannot add child to a Sequence (tip)')

    def __str__(self):
        return "Sequence {0.name!s} begins {1!s}...".format(self, self.seq[:10])


class Edge(object):
    """
    Class for the edges of the tree
    """

    def __init__(self, child, parent, length=None):
        self.child = child
        self.parent = parent
        self.length = length

    @property
    def child(self):
        return self._child

    @child.setter
    def child(self, new_child):
        self._child = new_child

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, new_parent):
        self._parent = new_parent
        new_parent.append_edge(self)

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, new_len):
        self._length = new_len

    @property
    def name(self):
        return 'Edge_'+ self.child.name + '_' + self.parent.name

    
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
                seqs.append(Sequence(seqid,seq, index=count))

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

    return dist/l


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

    genetic_dist = pd.DataFrame(columns=labels, data=dist_matrix, index=labels)

    genetic_dist.to_csv('hw3_genetic_distance.txt', sep='\t', float_format='%.15f')

    return genetic_dist

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

    dist_matrix = build_matrix([nodeslist[key] for key in nodeslist.keys()])

    count = nodeNum = len(nodeslist)
    
    while count>2:
        nodeNum += 1

        """ Calculate (joining) the Q-matrix """
        qmatrix = calc_qmat(dist_matrix)

        """ Branch length calculation """

        # Find the indices of the smallest (most negative) Q-value
        indx1,indx2 = np.unravel_index(qmatrix.values.argmin(), qmatrix.shape)
        lbl1 = qmatrix.columns[indx1]
        node1 = nodeslist[lbl1]
        lbl2 = qmatrix.columns[indx2]
        node2 = nodeslist[lbl2]

        # Calculate the distances to the new node
        deltaindx1 = (dist_matrix.iloc[indx1,indx2]/2 + 1/(2*(count-2)) *
                    (np.sum(dist_matrix.iloc[indx1,:]) - np.sum(dist_matrix.iloc[indx2,:])))

        deltaindx2 = dist_matrix.iloc[indx1,indx2] - deltaindx1
        
        assert deltaindx1>=0 and deltaindx2>=0, 'An edge length is less than zero'

        # Create and append edges and nodes
        newNode = 'Node'+str(nodeNum)
        nodeslist[newNode] = Node(name=newNode, index=nodeNum)
        newEdge1 = Edge(node1, nodeslist[newNode], length=deltaindx1)
        node1.parent = nodeslist[newNode]
        nodeslist[newNode].add_child(node1)
        edgelist[newEdge1.name] = newEdge1
        newEdge2 = Edge(node2, nodeslist[newNode], length=deltaindx2)
        node2.parent = nodeslist[newNode]
        nodeslist[newNode].add_child(node2)
        edgelist[newEdge2.name] = newEdge2

        """ Update distance matrix """
        new_node_frame = 1/2 * (dist_matrix[lbl1]+
                                dist_matrix[lbl2]-
                                dist_matrix[lbl1][lbl2])

        new_node_frame.name = newNode
        new_node_frame[newNode] = 0
        dist_matrix = dist_matrix.append(new_node_frame)
        dist_matrix[newNode] = new_node_frame

        dist_matrix.drop([lbl1, lbl2], axis=0, inplace=True)
        dist_matrix.drop([lbl1, lbl2], axis=1, inplace=True)

        count -= 1

    lbl1 = dist_matrix.columns[0]
    lbl2 = dist_matrix.columns[1]
    node1 = nodeslist[lbl1]
    node2 = nodeslist[lbl2]
    lastEdge = Edge(node1, node2, length=dist_matrix[lbl1][lbl2])
    edgelist[lastEdge.name] = lastEdge
    node1.parent = node2
    node2.add_child(node1)

    return nodeslist, edgelist

def print_edges(nodes):
    """
    Print the list of edges to a tab delimited file.
    """

    # Find the root node
    for key in nodes.keys():
        if nodes[key].parent is None:
            adam = nodes[key]

    def preorder(root,f):

        if root.edges:
            for edge in root.edges:
                f.write("{0}\t{1}\t{2}\n".format(edge.parent.index, edge.child.index, edge.length))
                preorder(edge.child,f)
 
    with open('hw3_edges.txt', 'w') as f:
        preorder(adam, f)        

def print_newick(nodes):
    """
    Print the NEWICK file with postorder format.
    """

    # Find the root node
    for key in nodes.keys():
        if nodes[key].parent is None:
            adam = nodes[key]

    def postorder(root):

        if root.edges:
            nodelist = []
            for edge in root.edges:
                name = postorder(edge.child)
                length = edge.length
                nodelist.append('{!s}:{:0.11f}'.format(name,length))
            return '('+','.join(nodelist)+')'
        else:
            return root.name

    with open('hw3_newick_tree.txt','w') as f:
        f.write(postorder(adam)+';')          


if __name__ == '__main__':

    assert len(sys.argv)==2,'Exactly one argument required.'

    fna_file = sys.argv[1]

    nodes, edges = build_tree(fna_file)
    print_edges(nodes)
    print_newick(nodes)