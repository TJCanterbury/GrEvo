#!/usr/bin/env python3

""" Align 2 networks where symmetry has been encoded in node alnmes """

__appalnme__ = 'Glob4.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.4'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from nxcode import readx
import numpy as np
import random

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

def add_ran_node(graph, reflect=False, middle=False):
	""" Add node """
	node = str(graph.number_of_nodes() + 1)
	while graph.has_node(node):
		node = str(int(node) + 1)
		
	if (np.random.randint(0,2) or reflect) and not middle:
		node = "L" + node
		graph.add_node(node)
		
		return graph, node

	else:
		graph.add_node(node)

		return graph, node

### Align ###
def S3(G1, G2, aln):
	edges_1 = set(G1.edges())
	edges_2 = set(G2.edges())
	c = 0
	d = 0
	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]):
			c += 1
	for f in edges_2:
		if G1.has_edge(f[0], f[1]):
			d += 1

	score = c / ( len(edges_1) - c + len(edges_2))

	return score

def Node_sim(G1, G2, node1, node2, aln):
	test_aln = aln.copy()
	test_aln[node1] = node2
	score = S3(G1, G2, aln)
	return score

def extend_aln(G1, G2, aln, n1, n2):
	""" attempt to add nodes to outer perimeter of subgraphs that 
	maintain isomorphism between these subgraphs and so building 
	the largest common subgraph """

	# for each of n1, align the node with the node in n2 with the highest similairty
	for node in n1:
		# build dictionary of similairites between node and all available 
		# nodes in the opposing graph
		if aln[node]:
			break
		Similarities = { node2: Node_sim(G1, G2, node, node2, aln) for node2 in n2 }
		if Similarities:
			# Assign alignment of first best node based on similarity scores
			aln[node] = max(Similarities, key=Similarities.get)
			aln[reflect_n(node)] = reflect_n(aln[node])

			# remove aligned nodes from subset that are unaligned
			n1.remove(node)
			try: # if node has a reflection it will be removed too
				n2.remove(aln[node])
			except:
				pass
			try: # if node has a reflection it will be removed too
				n1.remove(reflect_n(node))
			except:
				pass
			try: # if node has a reflection it will be removed too
				n2.remove(reflect_n(aln[node]))
			except:
				pass
	return aln

def ran_align(G1, G2, source, aln):
	""" randomly align nodes, with preference for those similar and within
	the same shortest path distance from the source/body """
	
	# make independent alignment from each node in source
	n1 = nx.single_source_shortest_path_length(G1, source[0][0])
	n2 = nx.single_source_shortest_path_length(G2, source[0][1])

	max_radius = min([max(n1.values()), max(n2.values())])
	radius = 1

	while radius <= max_radius:
		A_neighbs = [k for k in n1 if n1[k] == radius] 
		B_neighbs = [k for k in n2 if n2[k] == radius] 
		random.shuffle(A_neighbs)
		random.shuffle(B_neighbs)
		aln = extend_aln(G1, G2, aln, A_neighbs, B_neighbs)
		radius +=1
	
	n1 = [node for node in G1.nodes() if not aln[node]]
	n2 = [node for node in G2.nodes() if node not in aln.values()]
	aln = extend_aln(G1, G2, aln, n1, n2)

	# choose alignment with best symmetric substructure score
	return S3(G1, G2, aln), aln

def GrEvAl(Graph1, Graph2, source=[["body", "body"]]):
	""" Grow an alignment from source, 
	the known node correspondance such as the body """

	G1 = Graph1.copy()
	G2 = Graph2.copy()
	# Make the graphs the same size so that a global alignment can be made
	if G1.number_of_nodes() > G2.number_of_nodes():
		G1, G2 = G2, G1
	while G1.number_of_nodes() < G2.number_of_nodes():
		G1, node = add_ran_node(G1, middle=True)
	
	
	# start alignment with the known node correspondance(s) from source
	aln = dict.fromkeys(G1.nodes())
	for i in source:
		aln[i[0]] = i[1]
	
	score, aln = ran_align(G1.copy(), G2.copy(), source, aln)

	return score, aln

### Main ###
def main(argv):
	G1 = readx(argv[1])
	G2 = readx(argv[2])
	score, aln = Align(G1, G2, repeats=int(argv[3]), ret_aln=True)
	print(aln)
	print(score)
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)