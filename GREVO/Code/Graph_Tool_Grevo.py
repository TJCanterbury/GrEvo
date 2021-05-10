#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measure distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'GrEvo2.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
from numpy import genfromtxt
from time import time
import graph_tool.all as gt 
import graph_tool.draw as dr 
import graph_tool.generation as rg 
import graph_tool.util as ug 

## Functions ##
### General Graph Functions ###
def num_edges(mean = 2.5):
    """ Returns from a normal distribution an integer for the
    number of edges to add """
    randomInts = np.random.normal(loc=mean, size=1).astype(int)
    while randomInts <= 0:
        randomInts = np.random.normal(loc=mean, size=1).astype(int)
    return int(randomInts)

def correct_L_or_R(G):
    """ returns 0 for middle, 1 for left, 2 for right """
    for n in G.nodes():
        if n[-1:] == 'l':
            # If on left side add copy on right side
            n2 = 'L' + n[:-1]
            mapping = {n:n2}
            G = nx.relabel_nodes(G, mapping)
        
        # Check if on right side
        elif n[-1:] == 'r':
            # If on right side add copy on left side
            n2 = 'R' + n[:-1]
            print (n2)
            mapping = {n:n2}
            G = nx.relabel_nodes(G, mapping)
        
    return G

def prop_test(G):
    V = np.random.choice(list(G.vertices()))
    return V

def reflect_n(graph, node):
    """ returns the reflected node label """

    name = graph.vp.names[node]

    if name[0] == 'L':
        # If on left side add copy on right side
        name = 'R' + name[1:]
    
    elif name[0] == 'R':
        # If on right side add copy on left side
        name = 'L' + name[1:] 
    
    node = ug.find_vertex(graph, graph.vp.names, name)
    
    return node

def add_n_edges_to_node(graph, n, node):
    """ adds n random edges to a given node """
    target = graph.vertex(node).in_degree() + n
    c1 = list()
    c2 = list()
    if target > graph.num_vertices():
        return graph, c1, c2
    
    node2 = reflect_n(graph, node)
    print(node2)
    # Reflect all actions
    a1 = set(graph.get_all_neighbors(graph.vertex(node)))
    a2 = set(graph.get_all_neighbors(node2))
    
    total = set(list(graph.vertices()))

    while graph.get_total_degrees(node) < target:
        
        v = np.random.choice(total - a1)
        graph.add_edge(node, v)
        graph.add_edge(node2, reflect_n(graph, v))
    
    b1 = set(graph.get_all_neighbors(node))
    b2 = set(graph.get_all_neighbors(node2))
    
    c1 = list(b1 - a1)
    c2 = list(b2 - a2)
    
    return graph, c1, c2

def p_middle(nodes):
    """ gives proportion of nodes that are in the middle """
    total = len(nodes)
    x = 0
    for i in nodes:
        if i[-1:] == 'l':
            # If on left side add copy on right side
            x += 1
    return x/total

def make_middle(nodes):
    """ Returns True if character should be in middle """
    p = p_middle(nodes)
    roll = np.random.uniform()
    
    if roll <= p:
        return True
    
    return False

def add_ran_node(graph):
    """ Add node """
    al = list(graph.nodes())

    node = str(graph.number_of_nodes() + 1)

    if not make_middle(al):
        node = node + "l"
        graph.add_node(node)
        
        return graph, node

    else:
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
        node = np.random.choice(graph.get_vertices())
        #print(node)

    graph, b1, b2 = add_n_edges_to_node(graph, 1, node)

    #edges =  [node + ", " + ', '.join(b1), node2 + ", " +  ', '.join(b2)]
    #movement =  "char_grows- Added edges: " + edges[0] + "; " + edges[1]

    return graph #, movement

def char_shrinks(graph, node = None, with_removal = True):
    """plate shrinks, so loses edges. If leaf node remove """
    print(graph.nodes())
    if node == None:
        node = np.random.choice(graph.nodes())
    
    nodes = [node, reflect_n(node)]

    # if leaf node, nodes:
    if graph.degree(node) == 1 and with_removal:
        graph.remove_nodes_from(nodes)

        movement = "char_shrinks- node: " + ', '.join(nodes) + " removed"

        return graph, movement
    
    v = np.random.choice(list(graph.neighbors(node)))
    
    vs = [v, reflect_n(v)]
    edges = [(nodes[0], vs[0]), (nodes[1], vs[1])]

    # if leaf node, nodes:
    if graph.degree(node) == 1 and with_removal:
        graph.remove_nodes_from(nodes)

        movement = "char_shrinks- node: " + ', '.join(nodes) + " removed"

        return graph, movement
    

    else:
        # if not leaf node, remove all versions of edge from all versions of character:
        graph.remove_edges_from(edges)
        edges = [', '.join(list(elem)) for elem in edges]
        movement =  "char_shrinks- removed edges: " + edges[0] + "; " + edges[1]

        return graph, movement

def char_moves(graph, node = None):
    """ plate moves so edges are replaced """
    print(graph.nodes())
    if node == None:
        node = str(np.random.choice(graph.nodes()))

    movement = "char_moves- node: " + node
    
    graph, m1 = char_grows(graph, node)
    graph, m2 = char_shrinks(graph, node, with_removal=False)

    movement = "char_moves- node: " + node + \
        ". Edge change:" + m2 + " ---> " + m1

    return graph, movement

def char_gain(graph, Node = None):
    """ New plate emerges, new node with mean edges of 3 """
    print(graph.nodes())
    # Add node(s)
    graph, node1 = add_ran_node(graph)
    node2 = reflect_n(node1)
    graph.add_node(node2)
    
    if graph.number_of_nodes() > 2:
        # Add typical number of edges to node based on a normal pdf mean = 3
        n = num_edges()
        
        # reflected nodes
        graph, b1, b2 = add_n_edges_to_node(graph, n, node1)
        
        movement = "char_gain- node1: " + node1 + " new edges: " + ', '.join(b1) 
    
    else:
        graph.add_edge(node1, Node)
        
        movement = "char_gain- node1: " + node1 + " new edge: " + Node + "."
    
    return graph, movement

def char_loss(graph, Node = None):
    """ plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
    if Node == None:
        Node = str(np.random.choice(graph.nodes()))
    
    degree = graph.degree(Node)
    print(graph.nodes())
    if degree > 1:
        
        neighbs = np.random.choice(list(graph.neighbors(node)), num_edges())
        graph.add_edges_from(
            it.product(
                neighbs,
                neighbs
            )
        )
    
        graph.add_edges_from(
            it.product(
                reflect_n(neighbs),
                reflect_n(neighbs)
            )
        )
    
    if is_L_or_R(Node):
        graph.remove_node(reflect_n(Node))
    graph.remove_node(Node)
    
    movement = "char_loss- node: " + Node + ", " + reflect_n(Node) + "."

    return graph, movement

def char_merge(graph, u = None):
    """ 2 plates become one, 2 adjacent nodes become the same node, union of edges """
    print(graph.nodes())
    
    if u == None:
        u = str(np.random.choice(graph.nodes()))
    
    v = str(np.random.choice(list(graph.neighbors(u))))

    if not is_L_or_R(u) and not is_L_or_R(v): # 2 middle nodes merge
        
        print(1)
        print(u)
        print(reflect_n(u))
        print(v)
        print(reflect_n(v))
        graph = nx.contracted_nodes(graph, u, v)
        
        movement = "char_merge- Eater node: " + u + "." + " eaten node: " + v + "."
    
    elif not is_L_or_R(u) and is_L_or_R(v): # a middle node merges with 2 side nodes
        
        print(2)
        print(u)
        print(reflect_n(u))
        print(v)
        print(reflect_n(v))
        graph = nx.contracted_nodes(graph, u, v)
        graph = nx.contracted_nodes(graph, u, reflect_n(v))
        
        movement = "char_merge- Eater node: " + u + "." + " eaten nodes: " + v + ', ' +\
             reflect_n(v) + "."

    elif is_L_or_R(u) and not is_L_or_R(v): # 2 side nodes merge with a middle node
        
        print(3)
        print(u)
        print(reflect_n(u))
        print(v)
        print(reflect_n(v))
        graph = nx.contracted_nodes(graph, v, u)
        graph = nx.contracted_nodes(graph, v, reflect_n(u))
        
        movement = "char_merge- Eater node: " + v + ". " + \
            "eaten nodes: " + u + ', ' + reflect_n(u) + "."

    elif is_L_or_R(u) and is_L_or_R(v): # 2 side nodes merge, 
                                        # as do their relfections
        print(4)
        print(u)
        print(reflect_n(u))
        print(v)
        print(reflect_n(v))
        graph = nx.contracted_nodes(graph, u, v)
        graph = nx.contracted_nodes(graph, reflect_n(u), reflect_n(v))
    
        movement = \
            "char_merge- Eater nodes: " + u + ", " + reflect_n(u) + \
            ". Eaten nodes: " + v + ", " + reflect_n(v) + "."
    
    return graph, movement

def char_split(graph, u1 = None):
    """ 0ne plate becomes 2, half the instances of a node 
    are replaced with new node that will be adjacent to old node """
    print(graph.nodes())
    if u1 == None:
        u1 = str(np.random.choice(graph.nodes()))
    u2 = str(reflect_n(u1))
    graph, node = add_ran_node(graph)
    v1 = str(node)
    v2 = str(reflect_n(node))
    ud = int(graph.degree(u1) / 2)

    if ud > 1:
        # Collect 50% of adjacencies
        neighbs = set(np.random.choice(list(graph.neighbors(u1)), size=ud))
        neighbs = list(neighbs)
        for a in neighbs:

            a2 = reflect_n(a)
            # Give adjacencies to new node
            graph.add_edge(v1, a)
            graph.add_edge(v2, a2)

            # Remove adjacencies from old node0
            graph.remove_edge(u1, a)
            if u1 == a2:
                del a
                continue
            if is_L_or_R(u1):
                if is_L_or_R(a):
                    graph.remove_edge(u2, a2)
                elif not is_L_or_R(a):
                    graph.remove_edge(u2, a)

    # Add edge between split nodes
    graph.add_edge(u1, v1)
    graph.add_edge(u2, v2)
    
    movement = "char_split- from node(s): " + u1 + ', ' + u2 + " new node(s): " + v1 + ", " + v2 + "."

    return graph, movement

### Evo (hill climb) algorithm ###
def perturber(G, move = None, Node = None):
    """ Make random move """
    # Variables and Constraints
    G_Size = G.number_of_nodes()
    if Node == None:
        Node = str(np.random.choice(G.nodes()))
    if G_Size == 1:                                 # G can only gain or split
        move = np.random.choice([3, 6])
    if G_Size == G.degree(Node):                    # Node can't grow
        move = np.random.choice([1, 3, 4, 5, 6])
    if G.degree(Node) == 0:                         # Node must grow or be lost
        move = np.random.choice([0, 4, 5])
    if move == None:
        move = np.random.choice(range(6))
    
    # Move functions
    if move == 0:
        G, move = char_grows(G, Node)
    if move == 1:
        G, move = char_shrinks(G, Node)
    if move == 2:
        G, move = char_moves(G, Node)
    if move == 3:
        G, move = char_gain(G, Node)
    if move == 4:
        G, move = char_loss(G, Node)
    if move == 5:
        G, move = char_merge(G, Node)
    if move == 6:
        G, move = char_split(G, Node)

    return G, move 

def intelligence(G1, G2):
    """ Directs the perturber towards nodes and edges
    that are most wrong and towards the moves that are
    most likely to reduce the distance """

    return 0

def recorder(move):
    """ Store moves that improve portrait divergence score for future optimisation """
    return

def measurer(G1, G2):
    """ Measure Portrait divergeance (aiming for 0) """
    return

def climber(G1, G2, old_score = 1, start_time = None):
    """ Apply perturber, if score improved record and recurse """
    if start_time == None:
        end_time = time() + (5 * 60)
    morph, move = perturber(G1)
    score = 0.5#measurer(morph, G2)
    print(score)
    if score == 0:
        return morph
    elif end_time == time():
        #recorder(move)
        return morph

    elif score <= old_score:
        #recorder(move)
        climber(morph, G2, score, end_time)

### Business End ###
def main(argv):
        
    G1 = gt.Graph(directed=False) 
    e = np.genfromtxt(argv[1], dtype=str, delimiter=" ")
    G1.vp.names = G1.add_edge_list(e, hashed=True)
    name = G1.vp.names[5]
    name2 = reflect_n(G1, 5)
    #print(name)
    #print(name2)
    i = ug.find_vertex(G1, G1.vp.names, name2)
   # print(i)
    #print(G1.vertices())
    G1 = char_grows(G1)
    dr.graph_draw(G1, vertex_text=G1.vp.names)
    
    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)