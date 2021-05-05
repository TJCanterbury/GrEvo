#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measure distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'grevo.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import networkx as nx
from nxcode import readx
from nxcode import drawx

## Functions ##
### General Graph Functions ###
def num_edges(mean = 4):
    """ Returns from a normal distribution an integer for the
    number of edges to add """
    randomInts = np.random.normal(loc=mean, size=1).astype(int)
    while randomInts <= 0:
        randomInts = np.random.normal(loc=mean, size=1).astype(int)
    return int(randomInts)

def add_n_edges_to_node(graph, n, node):
    """ adds n random edges to a given node """
    target = graph.degree(node) + n
    a = set(graph.neighbors(node))
    while graph.degree(node) < target:
        v = np.random.choice(list(nx.non_neighbors(graph, node)))
        graph.add_edge(node, str(v))

    b = set(graph.neighbors(node))
    c = b - a
    return graph, c

def add_ran_node(graph):
    """ Add node """
    al = list(graph.nodes())
    numbers = [ int(x) for x in al ]
    node = str(np.max(numbers) + 1)
    graph.add_node(node)

    return graph, node

def Node_dif(a, b):
    """ Number of nodes difference """
    nns = len(a.number_of_nodes()) - len(b.number_of_nodes())

    return nns

### Morph Functions ###
def char_grows(graph, node = None):
    """ plate grows, so node gains edges """
    if node == None:
        node = np.random.choice(graph.nodes())

    graph, b = add_n_edges_to_node(graph, 1, node)

    movement = "char_grows- node: " + node + " new edge: " + ''.join(list(b))
    print(movement)

    return graph

def char_shrinks(graph, node = None):
    """plate shrinks, so loses edges. If leaf node remove """
    if node == None:
        node = np.random.choice(graph.nodes())
    
    if graph.degree(node) == 1:
        graph.remove_node(node)
        
        movement = "char_shrinks- node: " + node + " removed"
        print(movement)

        return graph

    v = np.random.choice(list(graph.neighbors(node)))
    graph.remove_edge(node, v)

    movement = "char_shrinks- node: " + node + " removed adjacency: " + v
    print(movement)

    return graph

def char_moves(graph):
    """ plate moves so edges are replaced """
    node = np.random.choice(graph.nodes())

    movement = "char_moves- node: " + node
    print(movement)
    
    char_grows(graph, node)
    char_shrinks(graph, node)

    return graph

def char_gain(graph):
    """ New plate emerges, new node with mean edges of 3 """
    # Add node
    graph, node = add_ran_node(graph)

    # Add typical number of edges to node based on a normal pdf mean = 3
    n = num_edges()
    graph, b = add_n_edges_to_node(graph, n, node)
    
    movement = "char_gain- node: " + node + " new edges: " + ', '.join((list(b))) + "."
    print(movement)

    return graph

def char_loss(graph):
    """ plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
    node = np.random.choice(graph.nodes())
    degree = graph.degree(node)
    if degree > 1:
        graph.add_edges_from(
            it.product(
                graph.neighbors(node),
                graph.neighbors(node)
            )
        )
    graph.remove_node(node)
    
    movement = "char_loss- node: " + node
    print(movement)

    return graph

def char_merge(graph):
    """ 2 plates become one, 2 adjacent nodes become the same node, union of edges """
    u = np.random.choice(graph.nodes())
    v = np.random.choice(list(graph.neighbors(u)))
    graph  = nx.contracted_nodes(graph, u, v, self_loops=False)
    
    movement = "char_merge- node1: " + u + " node2: " + v
    print(movement)

    return graph

def char_split(graph):
    """ 0ne plate becomes 2, half the instances of a node 
    are replaced with new node that will be adjacent to old node """
    u = np.random.choice(graph.nodes())
    graph, v = add_ran_node(graph)
    ud = int(graph.degree(u) / 2)
    
    if ud > 1:
        # Collect 50% of adjacencies
        neighbs = set(np.random.choice(list(graph.neighbors(u)), size=ud))
        
        # Give adjacencies to new node
        for a in neighbs:
            graph.add_edge(v, a)

        # Remove adjacencies from old node
        for a in neighbs:
            graph.remove_edge(u, a)

    # Add edge between split nodes
    graph.add_edge(u, v)
    
    movement = "char_split- from node: " + u + " new node: " + v
    print(movement)

    return graph

### Evo (hill climb) algorithm ###

### Business End ###
def main(argv):
    a = readx(argv[1])
    G = "Nodes before: " + ', '.join(list(a.nodes)) + "."
    print(G)
    a = char_grows(a)
    a = char_shrinks(a)
    a = char_moves(a)
    a = char_gain(a)
    a = char_loss(a)
    a = char_merge(a)
    a = char_split(a)
    G = "Nodes after: " + ', '.join(list(a.nodes)) + "."
    print(G)
    drawx(a, "../pond/ame.png")
    
if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)