#!/usr/bin/env python3

""" Align 2 networks where symmetry has been encoded in node names """

__appname__ = 'SymAln.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from networkx.convert import from_edgelist
from nxcode import readx
from nxcode import drawx
import numpy as np

### Align and Measure Similarity ###
def reflect_n(node):
    """ returns the reflected node label """

    if node[0] == 'L':
        # If on left side add copy on right side
        node = 'R' + node[1:]
    
    elif node[0] == 'R':
        # If on right side add copy on left side
        node = 'L' + node[1:] 
    
    return node

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

### Align and LCN2 Similarity ###
def score_alignment(G1, G2, aln):
	""" if subgraphs given by aln are isomorphic return true and the number of nodes aligned """
	A_nodes = [i[0] for i in aln]
	B_nodes = [i[1] for i in aln]

	A = G1.subgraph(A_nodes)
	B = G2.subgraph(B_nodes)
	A_nodes = A.number_of_nodes()
	if nx.is_isomorphic(A, B):
		return True, A_nodes
	else:
		return False, A_nodes

def grow_alignment(G1, G2, radius, aln = [["body", "body"]]):
	""" add random node to alignment """
	A_nodes = [i[0] for i in aln]
	B_nodes = [i[1] for i in aln]
	
	G1_nodes = nx.ego_graph(G1, "body", radius=radius, undirected= True).nodes()
	G2_nodes = nx.ego_graph(G2, "body", radius=radius, undirected= True).nodes()
	
	n_type = np.random.choice([0,1,2])
	A_neg_nod = [ node for node in G1_nodes if is_L_or_R(node) == n_type and node not in A_nodes]
	B_neg_nod = [ node for node in G2_nodes if is_L_or_R(node) == n_type and node not in B_nodes]
	
	if A_neg_nod and B_neg_nod:
		va = np.random.choice(A_neg_nod)
		vb = np.random.choice(B_neg_nod)
		aln.append([va, vb])
	
		if is_L_or_R(va):
			var = reflect_n(va)
			vbr = reflect_n(vb)
			aln.append([var, vbr])

	# Check whether the aligned nodes represent the same subgraph
	isomorphic, score = score_alignment(G1, G2, aln)
	if isomorphic:
		return aln, score
	return None, None

def common_neighbourhood(G1, G2, radius, node = "body"):
	A = nx.ego_graph(G1, node, radius=radius, undirected= True)
	B = nx.ego_graph(G2, node, radius=radius, undirected= True)
	if nx.is_isomorphic(A, B):
		return True, A.number_of_nodes()
	else:
		return False, A.number_of_nodes()

def LCN(G1, G2):
	Num = G2.number_of_nodes() 
	hmm = True
	radius = 0
	nodes = 1

	while hmm == True:
		old_nodes = nodes
		radius += 1
		hmm, nodes = common_neighbourhood(G1, G2, radius = radius)
	score = abs(1 - (old_nodes/Num))
	return score

def LCN2(G1, G2, aln = [["body", "body"]], stuck=0):
	Num = G2.number_of_nodes() 
	i = 1
	radius = 1
	score = None
	while i <= 1000 or score == None:
		possible_aln, Poss_score = grow_alignment(G1, G2, radius, aln.copy())
		if Poss_score:
			aln = possible_aln
			score = Poss_score
			radius += 1
		i+=1
	print(score)
	score = abs(1 - (score/Num))
	return aln, score
	
### Align and measurer Similarity ###
def LCN4(G1, G2, stuck=False, node="body"):
	""" Build an alignment and score it by the largest common subgraph """
	# Initiate variables
	NumA = G1.number_of_edges() 
	NumB = G2.number_of_edges()
	radius = 0
	neighbs_1 = nx.single_source_shortest_path_length(G1, "body")
	neighbs_2 = nx.single_source_shortest_path_length(G2, "body")
	max_radius_1 = max(neighbs_1.values())
	max_radius_2 = max(neighbs_2.values())

	if nx.is_isomorphic(G1, G2): # If solution found return score of 0
		return 0, radius

	while True: # Else, start by matching ego_graph about the body
		A = nx.ego_graph(G1, node, radius=radius, undirected= True)
		B = nx.ego_graph(G2, node, radius=radius, undirected= True)
		
		if nx.is_isomorphic(A, B):
			Asp = list(A.nodes())
			Bsp = list(B.nodes())
			radius += 1
		else:
			save_rad = (radius - 1)
			break
	
	while radius != max_radius_1:
		Asp, Bsp = extend_pos(G1, G2, radius, Asp, Bsp, neighbs_1, neighbs_2)
		radius += 1

	# Negative direction subgraph alignment
	rad1, rad2 = max_radius_1, max_radius_2
	Asn = []
	Bsn = []
	while rad1 >= 0 and rad2 >= 0:
		Asn, Bsn = extend_neg(G1, G2, rad1, rad2, Asn, Bsn, neighbs_1, neighbs_2)
		rad1 -= 1
		rad2 -= 1
	
	neg_size = len(set(Asp))
	pos_size = len(set(Asn))

	if neg_size >= pos_size: # which method gave largest common subgraph?
		Sub_n = neg_size
		Bs = Bsp
	else:
		Sub_n = pos_size
		Bs = Bsn
	
	# Score as proportion of most nodes aligned 
	# against the number of nodes in the biggest of the 2 graphs
	if NumB >= NumA: 
		score = abs(1 - (Sub_n / NumB))
	else:
		score = abs(1 - (Sub_n / NumA))

	return score, G2.subgraph(set(Bs)) #, save_rad

def extend_pos(G1, G2, radius, As, Bs, n1, n2):
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

def extend_neg(G1, G2, rad1, rad2, As, Bs, n1, n2):
	""" attempt to add nodes to outer perimeter of subgraphs that 
	maintain isomorphism between these subgraphs and so building 
	the largest common subgraph """

	A_neighbs = [k for k in n1 if n1[k] == rad1]
	B_neighbs = [k for k in n2 if n2[k] == rad2]
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

def main(argv):
	G1 = readx(argv[1])
	G2 = readx(argv[2])
	score, Sub_G1 = LCN4(G1, G2, stuck=True)
	drawx(Sub_G1, "display_test.png")
	print(score)
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)