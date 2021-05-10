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
from portrait_divergence import portrait_divergence

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
        
def is_L_or_R(node):
    """ returns 0 for middle, 1 for left, 2 for right """
    if node[0] == 'L':
        # If on left side add copy on right side
        return 1
    
    # Check if on right side
    elif node[0] == 'R':
        # If on right side add copy on left side
       return 2
    
    else:
        return 0

def reflect_n(node):
    """ returns the reflected node label """

    if node[0] == 'L':
        # If on left side add copy on right side
        node = 'R' + node[1:]
    
    elif node[0] == 'R':
        # If on right side add copy on left side
        node = 'L' + node[1:] 
    
    return node

def reflect_e(edgelist):
    """ Make edge list for symmetrical grevo derived graph """
    reflections = set()
    G_set = set( frozenset(element) for element in edgelist if element[0]!=element[1])
    # Iterate through nodes and add to list their reflections
    for edge in edgelist:
        reflection = []
        for node in edge:

            # Check if on left side
            if is_L_or_R(node) == 1:
                # If on left side add copy on right side
                reflection.append(node[:-1] + 'r')
            
            # Check if on right side
            elif is_L_or_R(node) == 2:
                # If on right side add copy on left side
                reflection.append(node[:-1] + 'l')

            # add any middle nodes
            else:
                reflection.append(node)
        if reflection[0] != reflection[1]:
            # Add each reflection to the reflections set
            reflections.add(frozenset(reflection))
    
    # Create set of symmetrical graph
    G_set = G_set.union(reflections)
    #convert to list for writing to file
    reflected_elist = [ list(element) for element in G_set ]

    return reflected_elist

def add_n_edges_to_node(graph, n, node):
    """ adds n random edges to a given node """
    target = graph.degree(node) + n
    c1 = list()
    c2 = list()
    if target > graph.number_of_nodes()-1:
        return graph, c1, c2
    node2 = reflect_n(node)
    
    # Reflect all actions
    a1 = set(graph.neighbors(node))
    a2 = set(graph.neighbors(node2))
    
    while graph.degree(node) < target:
        v = np.random.choice(list(nx.non_neighbors(graph, node)))
        graph.add_edge(node, str(v))
        graph.add_edge(node2, reflect_n(str(v)))
    
    b1 = set(graph.neighbors(node))
    b2 = set(graph.neighbors(node2))
    
    c1 = list(b1 - a1)
    c2 = list(b2 - a2)
    
    return graph, c1, c2

def add_ran_node(graph):
    """ Add node """
    al = list(graph.nodes())

    node = str(graph.number_of_nodes() + 1)

    if np.random.randint(0,2):
        node = "L" + node
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
        node = np.random.choice(graph.nodes())

    node2 = reflect_n(node)

    graph, b1, b2 = add_n_edges_to_node(graph, 1, node)

    edges =  [node + " - " + ', '.join(b1), node2 + " - " +  ', '.join(b2)]
    movement =  "char_grows- Added edges: " + edges[0] + "; " + edges[1]

    return graph, movement

def char_shrinks(graph, node = None):
    """plate shrinks, so loses edges. If leaf node remove """
    if node == None:
        node = np.random.choice(graph.nodes())
    
    # Pick neighbours of node, symmetrically
    nodes = [node, reflect_n(node)]
    v = np.random.choice(list(graph.neighbors(node)))
    vs = [v, reflect_n(v)]
    edges = [(nodes[0], vs[0]), (nodes[1], vs[1])]
    
    #remove collected edges
    graph.remove_edges_from(edges)
    edges = [' - '.join(list(elem)) for elem in edges]
    movement =  "char_shrinks- removed edges: " + edges[0] + "; " + edges[1]
    
    # If node no longer connected remove
    if graph.degree(node) == 0:
        graph.remove_nodes_from(nodes)

    return graph, movement

def char_moves(graph, node = None):
    """ plate moves so edges are replaced """
    if node == None:
        node = str(np.random.choice(graph.nodes()))
    
    graph, m1 = char_grows(graph, node)
    graph, m2 = char_shrinks(graph, node)

    movement = "char_moves- node: " + node + \
        ". Edge change:" + m2 + " ---> " + m1

    return graph, movement

def char_gain(graph, Node = None):
    """ New plate emerges, new node with mean edges of 3 """
    # Add node(s)
    graph, node1 = add_ran_node(graph)
    node2 = reflect_n(node1)
    graph.add_node(node1)
    graph.add_node(node2)
    graph.add_edge(node1, Node)
    graph.add_edge(node2, reflect_n(Node))
    
    if graph.number_of_nodes() > 2:
        # Add typical number of edges to node based on a normal pdf mean = 3
        n = num_edges()
        
        # reflected nodes
        graph, b1, b2 = add_n_edges_to_node(graph, n, node1)
        
        Node += ', '.join(b1)
    
    movement = "char_gain- node1: " + node1 + " new edge: " + Node + "."
    
    return graph, movement

def char_loss(graph, Node = None):
    """ plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
    if Node == None:
        Node = str(np.random.choice(graph.nodes()))
    
    degree = graph.degree(Node)

    if degree > 1:
        
        neighbs = np.random.choice(list(graph.neighbors(Node)), num_edges())
        ref_neighbs = [reflect_n(node) for node in neighbs]
        graph.add_edges_from(
            it.product(
                neighbs,
                neighbs
            )
        )
    
        graph.add_edges_from(
            it.product(
                ref_neighbs,
                ref_neighbs
            )
        )
    
    if is_L_or_R(Node):
        graph.remove_node(reflect_n(Node))
    graph.remove_node(Node)
    
    movement = "char_loss- node: " + Node + ", " + reflect_n(Node) + "."

    return graph, movement

def char_merge(graph, u = None):
    """ 2 plates become one, 2 adjacent nodes become the same node, union of edges """
    
    if u == None:
        u = str(np.random.choice(graph.nodes()))
    
    v = str(np.random.choice(list(graph.neighbors(u))))

    if u == reflect_n(v):
        node = u[1:]
        graph.add_node(node)
        graph = nx.contracted_nodes(graph, node, u)
        graph = nx.contracted_nodes(graph, node, v)
        
        movement = "char_merge- Eater node: " + node + "." + \
            " eaten node: " + u + ", " + v + "."
    

    if not is_L_or_R(u) and not is_L_or_R(v): # 2 middle nodes merge
        graph = nx.contracted_nodes(graph, u, v)
        
        movement = "char_merge- Eater node: " + u + "." + " eaten node: " + v + "."
    
    elif not is_L_or_R(u) and is_L_or_R(v): # a middle node merges with 2 side nodes
        graph = nx.contracted_nodes(graph, u, v)
        graph = nx.contracted_nodes(graph, u, reflect_n(v))
        
        movement = "char_merge- Eater node: " + u + "." + " eaten nodes: " + v + ', ' +\
             reflect_n(v) + "."

    elif is_L_or_R(u) and not is_L_or_R(v): # 2 side nodes merge with a middle node
        graph = nx.contracted_nodes(graph, v, u)
        graph = nx.contracted_nodes(graph, v, reflect_n(u))
        
        movement = "char_merge- Eater node: " + v + ". " + \
            "eaten nodes: " + u + ', ' + reflect_n(u) + "."

    elif is_L_or_R(u) and is_L_or_R(v): # 2 side nodes merge, 
                                        # as do their relfections
        graph = nx.contracted_nodes(graph, u, v)
        graph = nx.contracted_nodes(graph, reflect_n(u), reflect_n(v))
    
        movement = \
            "char_merge- Eater nodes: " + u + ", " + reflect_n(u) + \
            ". Eaten nodes: " + v + ", " + reflect_n(v) + "."
    
        
    return graph, movement

def char_split(graph, u1 = None):
    """ 0ne plate becomes 2, half the instances of a node 
    are replaced with new node that will be adjacent to old node """
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
            try:
                graph.remove_edge(u1, a)
                graph.remove_edge(u2, a2)
            except:
                pass

    # Add edge between split nodes
    graph.add_edge(u1, v1)
    graph.add_edge(u2, v2)
    
    movement = "char_split- from node(s): " + u1 + ', ' + u2 + " new node(s): " + v1 + ", " + v2 + "."

    return graph, movement

### Evo (hill climb) algorithm ###
def perturber(G1, move = None, Node = None):
    """ Make random move """
    #Deep copy of graph:
    G = G1.copy()

    #Variables and Constraints
    G_Size = G.number_of_nodes()
    
    if Node == None:
        Node = str(np.random.choice(G.nodes()))
    if Node[:4] == "body":
        return G, move 
    elif G_Size <= 2:              # G can only gain or split
        move = np.random.choice([3, 6])
    elif G_Size == G.degree(Node): # Node can't grow
        move = np.random.choice([1, 3, 4, 5, 6])
    elif G.degree(Node) == 0:      # Node must grow or be lost
        move = np.random.choice([0, 4, 5])
    elif move == None:
        move = np.random.choice(range(6))
    
    if move == 0: # grows
        G, move = char_grows(G, Node)
    if move == 1: # shrinks
        G, move = char_shrinks(G, Node)
    if move == 2: # moves
        G, move = char_moves(G, Node)
    if move == 3: # gain
        G, move = char_gain(G, Node)
    if move == 4: # loss
        G, move = char_loss(G, Node)
    if move == 5: # merge
        G, move = char_merge(G, Node)
    if move == 6: # split
        G, move = char_split(G, Node)

    try: # remove solitary nodes
        G.remove_edges_from(nx.selfloop_edges(G)) 
        solitary=[ n for n,d in G.degree_iter(with_labels=True) if d==0 ]
        G.delete_nodes_from(solitary)
    except:
        pass

    return G, move 

def intelligence(G1, G2):
    """ Directs the perturber towards nodes and edges
    that are most wrong and towards the moves that are
    most likely to reduce the distance """

    morph, move = perturber(G1)
    try:
        score = portrait_divergence(morph, G2)
    except:
        pass

    return 0

def measurer(G1, G2):
    eg1 = nx.convert_node_labels_to_integers(G1)
    eg2 = nx.convert_node_labels_to_integers(G2)
    nx.write_edgelist(eg1, "../nupond_res/eg1.txt", delimiter=" ")
    nx.write_edgelist(eg1, "../nupond_res/eg2.txt", delimiter=" ")

def recorder(move, score, print_ = True):
    """ Store moves that improve portrait divergence score for future optimisation """
    Line = move + " Score: " + str(score)+"\n"

    if print_ == True:
        print(Line)
    return Line

def climber(G1, G2, max_moves = 1000, sample_size = 100, file = "test.txt"):
    """ Apply perturber, if score improved record and recurse """
    Generation = 0
    score = 1
    moves = []
    graphs = [G1]
    old_score = portrait_divergence(G1, G2)
    move = "Start Distance"
    moves.append(recorder(move, old_score))
    size = sample_size
    bad_directions = []
    
    while Generation < max_moves:
        breadth = 0
        best_G = None
        best_move = None
        best_score = 1

        while breadth < size:
            try:
                morph, move = perturber(graphs[Generation])
                score = portrait_divergence(morph, G2)
                for i in bad_directions:
                    if nx.is_isomorphic(morph, i):
                        print("nope")
                        morph, move = perturber(graphs[Generation])
                        score = portrait_divergence(morph, G2)
                
                breadth += 1
            except: 
                pass   

            if score < best_score:
                best_score = score
                best_G = morph
                best_move = move

        if best_score < old_score:
            Generation += 1
            size = sample_size * Generation
            drawx(best_G, "display.png")
            graphs.append(best_G.copy())
            moves.append(recorder(best_move, best_score))
            old_score = best_score
        
            if best_score == 0:
                break

        else:
            
            size = sample_size
            Generation = 0
            score = 1
            graphs = [G1]
            old_score = portrait_divergence(G1, G2)
            move = "Start Distance"
            moves.append(recorder(move, old_score))
            bad_directions.append(morph.copy())

    
    return graphs, moves

### Business End ###
def main(argv):
    G1 = readx(argv[1])
    G2 = readx(argv[2])

    breadth = int(argv[3])
    Results_file = argv[4]
    
    G3, moves = climber(G1, G2, sample_size = breadth, file = Results_file)

    Gen= 0
    for i in G3:
        printer = "../nupond/" + "morphling_" + file[:-3] + \
              "_" + str(Gen) + ".png"
        drawx(i, printer)
        Gen += 1

    for i in moves:
        with open(Results_file, "a") as myfile:
            myfile.write(i)

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)