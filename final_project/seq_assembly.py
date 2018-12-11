###
# CSCI 5481 Final Project
#
# Prof. Knights
#
# Author: Christopher White
# Date: December 11, 2018
###

import Bio
from Bio import SeqIO
from random import randrange, choice
from collections import defaultdict
from copy import deepcopy
from graphviz import Digraph
from difflib import SequenceMatcher


class DirectedGraph(object):
    """
    This idea implements an directed graph. The
    vertices and edges are added with instance methods. The graph
    can be plotted as well.
    """

    def __init__(self, k=1):
        
        self.vertices = defaultdict(int)
        self.edges = defaultdict(list)
        self.k = k

    def add_vertex(self, vertex):
        """
        Adds the vertex and increments the count on the vertex
        """
        if vertex not in self.vertices.keys():
            self.vertices[vertex] = 0

    def add_edge(self, edge):
        """
        Adds the edge to the dictionary of edges.  An *edge*
        is a tuple that has a *tail* and *head*.
        """
        tail, head = edge

        self.vertices[tail] -= 1
        self.vertices[head] += 1
        self.edges[tail].append(head)

    @staticmethod
    def count_edges(edges):
        """
        This function counts the number of edges in the
        edges dict *edges*
        """

        return sum([len(edge_list) for edge_list in edges.values()])

    def count_vertices(self, vertex, edges, visited):
        """
        This function counts the reachable vertices from a
        *vertex*. This is used in the bridge finding.
        Taken from:
        https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/
        """

        count = 1
        visited[vertex] = True

        for v in edges[vertex]:
            if visited[v] == False:
                count += self.count_vertices(v, edges, visited)

        return count
     
    def plot_graph(self, name='Graph'):
        graph = Digraph(name, format='png')
        for (tail, heads) in self.edges.items():
            for head in heads:
                graph.edge(tail, head)

        graph.view()

class EulerianGraph(DirectedGraph):
    """
    Some of these ideas were found online though the code is
    original. The website I reference is:
    https://www.geeksforgeeks.org/euler-circuit-directed-graph/
    """

    def __init__(self, base_list=None, k=1):
        super().__init__(k=k)

        self.path = []

        if base_list:
            self.add_sequence(base_list)

    def add_sequence(self, base_list):
        """
        Iterate over the sequence in *base_list* create the k-mers
        and adding them to the edges and vertices.
        """

        k = self.k
        L = len(base_list)

        for n in range(L-k+1):
            vertex = base_list[n:n+k]

            if n+k<L:
                tail, head = vertex, base_list[n+1:n+k+1]
                self.add_edge((tail, head))
            else:
                self.add_vertex(vertex)

    def valid_edge(self, edge, edges, vertices):
        """
        Determines if an edge is valid to traverse in the Euler path.
        Taken from:
        https://www.geeksforgeeks.org/fleurys-algorithm-for-printing-eulerian-path/
        """
        tail, head = edge

        if len(edges[tail]) == 1:
            # If this is the only edge, then it is valid.
            return True
        else:
            # Otherwise, check if the edge is a bridge.
            visited = dict(zip(vertices.keys(), [False]*len(vertices)))
            
            max_vertices = self.count_vertices(head, edges, visited)
            
            edges[tail].remove(head)
            visited = dict(zip(vertices.keys(), [False]*len(vertices)))
            num_vertices = self.count_vertices(head, edges, visited)
            edges[tail].append(head)

            if max_vertices > num_vertices:
                return False
            else:
                return True

    def find_path(self):
        self.path=[]
        start_value = min(self.vertices.values())

        start_vertices = [key for (key,value) in self.vertices.items() if value == start_value]

        edges = deepcopy(self.edges)
        vertices = deepcopy(self.vertices)
        
        start = start_vertices[0]
        max_verts = 0
        for vertex in start_vertices:
            # Find the most reachable vertices from each possible start
            visited = dict(zip(self.vertices.keys(), [False]*len(self.vertices)))
            count = self.count_vertices(vertex, edges, visited)

            if count > max_verts:
                max_verts = count
                start = vertex

        self.path = []

        self._euler(start, edges, vertices)

        return self.path

    def _euler(self, start, edges, vertices):
        """
        A recursive function that tries each edge and checks if 
        the all edges are used. If they are all used, the path
        is generated. Otherwise, no path is found.
        """

        for node in edges[start]:
            edge = (start, node)
            if self.valid_edge(edge, edges, vertices):
                edges[start].remove(node)
                if self.count_edges(edges) > 0:
                    if self._euler(node, edges, vertices):
                        # The final node was reached so append
                        # the edge to the path
                        self.path.insert(0, edge)
                        return True
                    else:
                        edges[start].append(node)
                else:
                    # The final edge was used so the path is complete
                    self.path.append(edge)
                    return True

        return False

            

    def assemble(self):
        sequence, _ = self.path[0]
        for _, head in self.path:
            sequence += head[-1]

        self.sequence = sequence


class HamiltionianGraph(DirectedGraph):
    """
    Hamiltonian graph of reads using the basic directed graph
    as the basis for the class. A minimum overlap is used to
    break a sequence into k-mers.
    """

    def __init__(self, base_list=None, k=2, min_olap=1):
        super().__init__(k=k)
        self.min_olap = min_olap

        if base_list:
            self.add_sequence(base_list)

    def add_sequence(self, base_list):
        """
        Iterate over the sequence in *base_list* create the k-mers
        and adding them to the edges and vertices.
        """

        k = self.k
        L = len(base_list)

        for n in range(L-k+1):
            vertex = base_list[n:n+k]

            if n+k<L:
                tail, head = vertex, base_list[n+1:n+k+1]
                self.add_edge((tail, head))
            else:
                self.add_vertex(vertex)

        if L<k:
            self.add_vertex(vertex)

    def compute_overlap(self):
        """
        Take the edges dict and compute the overlap between
        every pair of vertices
        """

        overlap = defaultdict(dict)

        vertices = list(self.vertices.keys())

        for vertex in vertices:
            for node in vertices:
                if node==vertex:
                    continue
                else:
                    n = self._lcs(vertex, node)
                    m = 0
                    if n>0:
                        if node in overlap[vertex]:
                            m=overlap[vertex][node]
                        overlap[vertex][node] = max([m,n])
                    if n<0:
                        if vertex in overlap[node]:
                            m=overlap[node][vertex]
                        n *= -1
                        overlap[node][vertex] = max([m,n])
            
            vertices.remove(vertex)

        self.overlap = overlap

    def _lcs(self,a,b):
        """
        Function to compute the longest common substring.
        """

        s = SequenceMatcher(None)
        s.set_seq2(b)
        len_b = len(b)
        s.set_seq1(a)
        len_a = len(a)

        olaps = s.get_matching_blocks()

        if s.ratio() < self.min_olap/max([len_a,len_b]):
            return 0

        for (i,j,n) in olaps:
            if j==0 and i+n==len_a:
                # Then a is the tail and be is the head
                return n

            if i==0 and j+n==len_b:
                # Then b is the tail and a is the head
                return -1*n
     
        return 0

    def plot_graph(self, name='Graph'):
        graph = Digraph(name, format='png')
        for tail in self.overlap:
            for head, olap in self.overlap[tail].items():
                graph.edge(tail, head, label=str(olap))

        graph.view()

def parse_sequence(sequence):
    """
    Take in a SeqRecord and return a base sequence string 
    that has the gaps removed and all of the characters 
    changed to uppercase.
    """

    return str(sequence.seq).upper().replace('N','')


def pellets(sequence, length=100, errors=False):
    """
    Take a *sequence* and return a random section of *length* for
    each call to this function.

    Errors are added at about 1 per 100 bp if errors are enabled.
    """

    bases = 'AGCT'

    l = len(sequence)
    while True:
        idx = randrange(l-length+1)
        pellet = sequence[idx:idx+length]
        if errors:
            err_idx = randrange(100)
            if err_idx<length:
                err_char = pellet[err_idx]
                pellet[err_idx] = choice(bases.replace(err_char,''))

        yield pellet

def calc_match(seq1, seq2):
    """
    Take two sequences and calculate how closely they match.
    """

    matcher = SequenceMatcher(seq1, seq2)

    matcher.get_matching_blocks()

    return matcher.ratio()


if __name__=="__main__":

    sequence = SeqIO.read('chr1.fa','fasta')
    
    base_sequence = parse_sequence(sequence)

    pellet = pellets(base_sequence, length=25)

    seq_reads = [ next(pellet) for _ in range(len(base_sequence)//10) ]

    
