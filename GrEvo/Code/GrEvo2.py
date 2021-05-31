#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measurer distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'GrEvo2.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.2'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import networkx as nx
import gressure as gr
from nxcode import readx
from nxcode import drawx

## Classes ##
class Morphling(nx.Graph):
    def __init__(graph, data=None, **attr):
        nx.Graph.__init__(graph, data, **attr)
    
    ### Basic Morph Functions ### 
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
    
        v = "body"
        while v == 'body':
            v = np.random.choice(list(graph.neighbors(node)))
        vs = [v, reflect_n(v)]
        edges = [(nodes[0], vs[0]), (nodes[1], vs[1])]
        
        #remove collected edges
        graph.remove_edges_from(edges)
        edges = [' - '.join(list(elem)) for elem in edges]
        movement =  "char_shrinks- removed edges: " + edges[0] + "; " + edges[1]
    
        return graph, movement
    
    def char_moves(graph, node = None):
        """ plate moves so edges are replaced """
        if node == None:
            node = str(np.random.choice(graph.nodes()))
        
        graph, m1 = graph.char_grows(node)
        graph, m2 = graph.char_shrinks(node)
    
        movement = "char_moves- node: " + node + \
            ". Edge change: " + m2 + " ---> " + m1
    
        return graph, movement
    
    def char_gain(graph, Node = None):
        """ New plate emerges, new node with mean edges of 3 """
        graph, node1 = add_ran_node(graph)
        node2 = reflect_n(node1)
        graph.add_node(node1)
        graph.add_node(node2)
        graph.add_edge(node1, Node)
        graph.add_edge(node2, reflect_n(Node))
        
        if graph.number_of_nodes() > 2:
            n = graph.number_of_nodes()
            while n >= graph.number_of_nodes():
                # Add typical number of edges
                n = num_edges(graph)
            
            # reflected nodes
            graph, b1, b2 = add_n_edges_to_node(graph, n, node1)
            
            Node += ", " + ', '.join(b1) + ". Node_2: " + node2 +\
                 " new edge(s): " + reflect_n(Node) + ", " + ', '.join(b2)
        
        movement = "char_gain- Node_1: " + node1 + " new edge(s): " + Node + "."
        
        return graph, movement
    
    def char_loss(graph, Node = None):
        """ plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
        
        neighbs = nx.single_source_shortest_path_length(graph, 'body')
        neighbs = [k for k in neighbs if neighbs[k] > 1]
        if neighbs:
            Node = np.random.choice(neighbs)
        else:
            Node = np.random.choice(graph.nodes())
        degree = graph.degree(Node)
    
        if degree > 1:
            
            neighbs = np.random.choice(list(graph.neighbors(Node)), \
                np.random.randint(0, graph.number_of_nodes()))
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
        
        v = 'body'
        while v == 'body':
            v = str(np.random.choice(list(graph.neighbors(u))))
        
        if u == reflect_n(v):
            node = u[1:]
            graph.add_node(node)
            graph = nx.contracted_nodes(graph, node, u)
            graph = nx.contracted_nodes(graph, node, v)
            
            movement = "char_merge- Eater node: " + node + "." + \
                " eaten node: " + u + ", " + v + "."
    
        elif not is_L_or_R(u) and not is_L_or_R(v): # 2 middle nodes merge
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
            u1 = np.random.choice(graph.nodes())
        
        u2 = reflect_n(u1)
        graph, node = add_ran_node(graph)
        v1 = node
        v2 = reflect_n(node)
        ud = int(graph.degree(str(u1)))
    
        if ud > 1:
            # Collect 50% of adjacencies
            neighbs = set(np.random.choice(list(graph.neighbors(u1)), size=np.random.randint(1, ud)))
            neighbs = list(neighbs)
            for a in neighbs:
                a2 = reflect_n(a)
                # Give adjacencies to new node
                graph.add_edge(v1, a)
                graph.add_edge(v2, a2)
    
                # Remove adjacencies from old node0
                try:
                    graph.remove_edge(u1, a)
                except:
                    pass
                try:
                    graph.remove_edge(u2, a2)
                except:
                    pass
                
        if is_L_or_R(v1) and not is_L_or_R(u1):
            graph.add_edge(v1, v2)
            graph.remove_node(u1)
        else:
            # Add edge between split nodes
            graph.add_edge(u1, v1)
            graph.add_edge(u2, v2)
        
        movement = "char_split- from node(s): " + u1 + ', ' + u2 + " new node(s): " + v1 + ", " + v2 + "."
    
        return graph, movement
    
    def char_expansion(graph, u1 = None):
        """ Add node to edge """
        u2 = reflect_n(u1)
        # choose neighbour:
        n1 = np.random.choice(list(graph.neighbors(u1)), size=1)[0]
        n2 = reflect_n(n1)
        
        if is_L_or_R(n1) or is_L_or_R(u1):
            # New nodes:
            graph, v1 = add_ran_node(graph, reflect=True)
            v2 = reflect_n(v1)
            graph.add_node(v2)
    
            graph.add_edge(u1, v1)
            graph.add_edge(n1, v1)
            graph.add_edge(u2, v2)
            graph.add_edge(n2, v2)
            movement = "char_expansion- from edge(s): " + u1 + "-" + n1 + ', ' + u2 +  "-" + n2 + \
                " new node(s): " + v1 + ", " + v2 + "."
        
        else:
            # New nodes:
            graph, v1 = add_ran_node(graph, middle=True)
    
            graph.add_edge(u1, v1)
            graph.add_edge(n1, v1)
    
            movement = "char_expansion- from edge(s): " + u1 + "-" + n1 + ', ' + u2 +  "-" + n2 + \
                " new node(s): " + v1 + "."
    
    
        # Remove old edges
        try:
            graph.remove_edge(u1, n1)
        except:
            pass
        try:
            graph.remove_edge(u2, n2)
        except:
            pass
        
        
        return graph, movement
    
    def char_squeein(graph, node=None): #Incomplete
        return None
    
    def char_squeeout(graph, node=None): #Incomplete
        return None
    
    def seg_duplicate(graph, node=None): #Incomplete

        return 0

    def __str__(graph, file="display.png"):
        drawx(graph, file)
        return 0

## Functions ##
### General Graph Functions #####
def num_edges(G):
    """ Returns from a normal distribution an integer for the
    number of edges to add """
    my_degrees = G.degree()
    degree_values = [v for k, v in my_degrees]
    sum_G = sum(degree_values)

    mean = (sum_G / G.number_of_nodes()) - 1 # minus one because this will be used after one edge has already been added!
    randomInts = np.random.normal(loc=mean, size=1).astype(int)
    
    while randomInts < 0:
        randomInts = np.random.normal(loc=mean, size=1).astype(int)
    
    return int(randomInts)

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

    node = gr.Reflect([node])[0]
    
    return node

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
        # identify second degree neighbours to connect with
        neighbs = nx.single_source_shortest_path_length(graph, node, cutoff=2)
        second_neighbs = [k for k in neighbs if neighbs[k] == 2]

        if not second_neighbs:
            v = choose_ran_node(graph, exclude = [node, reflect_n(node)])
        else:
            # pick one
            v = np.random.choice(second_neighbs)

        # add edge between target node and 2nd degree neighbour, and reflect
        graph.add_edge(node, str(v))
        graph.add_edge(node2, reflect_n(str(v)))
    
    b1 = set(graph.neighbors(node))
    b2 = set(graph.neighbors(node2))
    
    c1 = list(b1 - a1)
    c2 = list(b2 - a2)
    
    return graph, c1, c2

def add_ran_node(graph, reflect=False, middle=False):
    """ Add node """
    node = str((graph.number_of_nodes() * 2) + 1)

    if (np.random.randint(0,2) or reflect) and not middle:
        node = "L" + node
        graph.add_node(node)
        
        return graph, node

    else:
        graph.add_node(node)

        return graph, node

def choose_ran_node(graph, exclude = ["body"]):
    """ choose random node that is not of a list of excluded nodes """
    v = exclude[0]

    while v in exclude:
        v = np.random.choice(graph.nodes())

    return v

### Align and measurer Similarity ###
def LCN4(G1, G2, stuck=False, node="body"):
    """ Build an alignment and score it by the largest common subgraph """
    # Initiate variables
    NumA = G1.number_of_nodes() 
    NumB = G2.number_of_nodes()
    radius = 0
    neighbs_1 = nx.single_source_shortest_path_length(G1, "body")
    neighbs_2 = nx.single_source_shortest_path_length(G2, "body")
    max_radius_1 = max(neighbs_1.values())
    
    if nx.is_isomorphic(G1, G2): # If solution found return score of 0
        return 0, radius

    while radius <= max_radius_1: # Else, start by matching ego_graph about the body
        A = nx.ego_graph(G1, node, radius=radius, undirected= True)
        B = nx.ego_graph(G2, node, radius=radius, undirected= True)
        
        if nx.is_isomorphic(A, B):
            As = list(A.nodes())
            Bs = list(B.nodes())
            save_rad = radius
            radius += 1
        else:
            break
    
    # extend the common connected subgraph as much as possible 
    done = False
    while not done:
        old = As.copy()
        As, Bs = extend_light(G1, G2, radius, As, Bs, neighbs_1, neighbs_2)
        radius += 1
    
        if old == As:
            done = True

    # Score as proportion of nodes aligned 
    # against the number of nodes in the biggest of the 2 graphs
    Sub_n = len(set(As))
    if NumB >= NumA: 
        score = abs(1 - (Sub_n / NumB))
    else:
        score = abs(1 - (Sub_n / NumA))

    return score, save_rad

def extend_light(G1, G2, radius, As, Bs, n1, n2):
    """ attempt to add nodes to outer perimeter of subgraphs that 
    maintain isomorphism between these subgraphs and so building 
    the largest common subgraph """

    A_neighbs = [k for k in n1 if n1[k] == radius]
    B_neighbs = [k for k in n2 if n2[k] == radius]
    old_As = As.copy()
    old_Bs = Bs.copy()

    for node_A in set(A_neighbs):
        for node_B in set(B_neighbs):
                
            old_As.append(node_A)
            old_As.append(reflect_n(node_A))
            old_Bs.append(node_B)
            old_Bs.append(reflect_n(node_B))
            
            A0 = G1.subgraph( set(old_As) )
            B0 = G2.subgraph( set(old_Bs) )

            if nx.is_isomorphic(A0, B0):
                As += node_A + reflect_n(node_A)
                Bs += node_B + reflect_n(node_B)
        
    return As, Bs

def measurer(G1, G2, stuck=False, use_rad=False):
    score, rad = LCN4(G1, G2, stuck=stuck)
    if use_rad:
        return score, rad
    return score

### Evo (hill climb) algorithm ###
def perturber(G1, move = None, Node = None, radius=1):
    """ Make random move """

    try_again = True
    while try_again == True:
        #Copy of graph:
        G = G1.copy()

        #Variables and Constraints
        G_Size = G.number_of_nodes()

        if Node == None:
            neighbs = nx.single_source_shortest_path_length(G, 'body')
            if radius == 0:
                radius += 1
            neighbs = [k for k in neighbs if neighbs[k] >= radius]
            Node = np.random.choice(neighbs)

        if move == None:
            move = np.random.choice(range(0, 8))
            if G_Size <= 2:              # G can only gain or split
                move = np.random.choice([3, 6, 7])
            elif G_Size == G.degree(Node): # Node can't grow
                move = np.random.choice([1, 3, 4, 5, 6, 7])
            elif G.degree(Node) == 0:      # Node must grow or be lost
                move = np.random.choice([0, 4, 5])

        if move == 0: # grows
            G, move = G.char_grows(Node)
        if move == 1: # shrinks
            G, move = G.char_shrinks(Node)
        if move == 2: # moves
            G, move = G.char_moves(Node)
        if move == 3: # gain
            G, move = G.char_gain(Node)
        if move == 4: # loss
            G, move = G.char_loss(Node)
        if move == 5: # merge
            G, move = G.char_merge(Node)
        if move == 6: # split
            G, move = G.char_split(Node)
        if move == 7: # Graph expansion (replace edge with node)
            G, move = G.char_expansion(Node)

        neighbs = nx.single_source_shortest_path_length(G, 'body')
        neighbs = [k for k in neighbs]
        
        # Ensure still attached to body
        if len(neighbs) > 1:
            try_again = False
        else:
            G = G1.copy()
            Node == None

    # remove selfloops
    G.remove_edges_from(nx.selfloop_edges(G))   

    # remove solitary nodes
    remove = [node for node,degree in dict(G.degree()).items() if degree == 0]
    G.remove_nodes_from(remove)
    
    return G, move 

def recorder(move, score, print_ = True):
    """ Store moves that improve portrait divergence score for future optimisation """
    Line = move + " Score: " + str(score)+"\n"

    if print_ == True:
        print(Line)
    return Line

def searcher(G1, G2, size, attempts, dead_ends, Generation, rad = 1):
    old_score = measurer(G1, G2)
    breadth = 0
    best_G = 0
    best_move = None
    best_score = float(1)
    stuck = 0
    trying = 0
    best_rad = rad

    while breadth < size or best_score >= old_score:
        
        # Make a random move and measurer the effect
        morph, move = perturber(G1, radius=best_rad)
        score, rad = measurer(G1=morph, G2=G2, stuck=stuck, use_rad=True)
        breadth += 1

        # Escape dead end and record its graph
        if stuck >= attempts:    
            return old_score, None, G1, rad
        
        # Check if graph is isomorphic with a previous dead end 
        # and if so ignore this attempt
        try_again = False 
        for i in dead_ends:
            if nx.is_isomorphic(morph, i):
                try_again = True  
        if try_again:
            trying += 1
            if trying < size:
                breadth -= 1
            continue
        
        # keep the latest best score/move/graph
        if score < best_score and not move == None:# and np.random.randint(0, Generation+2):
            best_score = score
            best_G = morph
            best_move = move

        if rad > best_rad:
            best_rad = rad
        
        
        # Record how stuck we are on this generation
        if breadth > size:
            stuck += 1
            breadth = 0
            print(stuck)
    
    return best_score, best_move, best_G, best_rad

def climber(G1, G2, sample_size = 100, file = "test.txt", attempts = 10, goal = 4):
    """ Apply perturber, if score improved record and recurse """
    Generation = 0
    moves1 = []
    g1s = [G1]
    old_score = measurer(G1, G2)
    move = "Start Distance"
    moves1.append("G1 " + recorder(move, old_score))
    dead_end1 = []
    size = sample_size
    best_score = old_score
    best_rad=1
    
    while best_score != 0:
        best_score, best_move, best_G, best_rad = searcher(G1=g1s[Generation], 
            G2=G2, size=size, attempts=attempts, dead_ends=dead_end1, 
            Generation=Generation, rad=best_rad)
        
        if best_score == 0:
            size = sample_size * (Generation + 1)
            best_G.__str__()
            g1s.append(best_G.copy())
            moves1.append(recorder(best_move, best_score))
            return g1s, moves1

        elif best_move and (Generation + 1) < goal:
            Generation += 1
            size = sample_size * (Generation + 1)
            best_G.__str__()
            g1s.append(best_G.copy())
            moves1.append(recorder(best_move, best_score))
        
        else:
            dead_end1.append(g1s[Generation]) 
            Generation = 0
            g1s = [G1]
            move = "Start Distance"
            moves1.append("G1 " + recorder(move, old_score))
            size = sample_size
            best_rad=1 

### Business End ###
def main(argv):
    G1 = readx(argv[1])
    G2 = readx(argv[2])
    
    G1 = Morphling(G1)
    
    if argv[3] == "manual":
        if len(argv) == 6:
            move=int(argv[5])
        else:
            move = None
        morph, move = perturber(G1, Node=argv[4], move=move)
        score = measurer(morph, G2)
        recorder(move, score)

        morph.__str__("display_test.png")
    
    else:
        breadth = int(argv[3])
        Results_file = argv[4]
        attempts = int(argv[5])
        goal = int(argv[6])

        G3, moves = climber(G1, G2, sample_size=breadth, file=Results_file, attempts=attempts, goal=goal)

        Gen= 0
        for i in G3:
            printer = "../Results/" + "morphling_" + Results_file + \
                  "_" + str(Gen) + ".png"
            drawx(i, printer)
            Gen += 1
        
        moves_results = "../Results/" + Results_file + "_moves.txt"
        for i in moves:
            with open(moves_results, "a") as myfile:
                myfile.write(i)
    
    return 0

if __name__ == "__main__": 
    """Makes sure the "main" function is called from command line"""  
    status = main(sys.argv)
    sys.exit(status)