#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measure distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'grev.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import networkx as nx
from nxcode import readx
from nxcode import drawx
import csv

## Functions ##
### General Graph Functions ###
def num_edges(mean = 2.5):
    """ Returns from a normal distribution an integer for the
    number of edges to add """
    randomInts = np.random.normal(loc=mean, size=1).astype(int)
    while randomInts <= 0:
        randomInts = np.random.normal(loc=mean, size=1).astype(int)
    return int(randomInts)

def reflect_e(edgelist):
    """ If nodes are left or right return full symmetrical set of nodes """
    reflections = {}
    G_set = set( frozenset(element) for element in G_list )
    # Iterate through nodes and add to list their reflections
    for edges in edgelist:
        for nodes in edges:
            reflection = []

            # Check if on left side
            if i[-1:] == l:
                # If on left side add copy on right side
                reflection.append(i[:-1] + r)
            
            # Check if on right side
            elif i[-1:] == r:
                # If on right side add copy on left side
                reflection.append(i[:-1] + l)
            
        reflections.add = frozenset(reflection)
    
    G_set = G_set.union(reflections)

    reflected_elist = [ list(element) for element in G_set ]

    return reflected_elist

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
    numbers = [ int(x) for x in al ] # wont work
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
    move = int(argv[2])
    if move == 1:
        a = char_grows(a)
    if move == 2:
        a = char_shrinks(a)
    if move == 4:
        a = char_moves(a)
    if move == 8:
        a = char_gain(a)
    if move == 16:
        a = char_loss(a)
    if move == 32:
        a = char_merge(a)
    if move == 64:
        a = char_split(a)
    G = "Nodes after: " + ', '.join(list(a.nodes)) + "."
    print(G)

    e_list = reflect_e(G.edgelist())
    with open("Changeling.txt", "wb") as f:
        writer = csv.writer(f, delimiter=' ')
        writer.writerows(e_list)
    
if __name__ == "__main__": 