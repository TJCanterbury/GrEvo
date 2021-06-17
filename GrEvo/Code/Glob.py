#!/usr/bin/env python3

""" Align 2 networks where symmetry has been encoded in node alnmes """

__appalnme__ = 'Glob.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from networkx.classes import graph
from Morphlings import Plates
import random
import pprint

### Align ###
def S3(G1, G2, aln):
	""" Calculate symmetric substructure score described by Saraph and Milenković
	(Saraph V. Milenković T. (2014) MAGNA: maximizing accuracy in global network alignment. Bioinformatics, 30, 2931–2940.) """
	edges_1 = set(G1.edges())
	edges_2 = set(G2.edges())
	c = 0
	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]) or G2.has_edge(aln[e[1]], aln[e[0]]):
			c += 1

	score = c / ( len(edges_1) - c + len(edges_2))
	
	return score

def ics_score(G1, G2, aln):
	""" Calculate Edge Correctness score for incomplete data """
	edges_1 = set(G1.edges())
	
	c = 0
	subnodes = []
	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]) or G2.has_edge(aln[e[1]], aln[e[0]]) \
			and G1.is_L_or_R(e[0]) == G1.is_L_or_R(aln(e[0])) \
			and G1.is_L_or_R(e[1]) == G1.is_L_or_R(aln(e[1])):
			c += 1
			subnodes.append(aln[e[0]])
			subnodes.append(aln[e[1]])
	subG1 = G2.subgraph(subnodes)
	sublen = subG1.size()
	score = c /( sublen + len(edges_1) - c)
	return score

def Node_sim(G1, G2, node1, node2, aln):
	""" Determine node similairty based on symmetry,
	and alignment """
	D_1 = G1.degree(node1)
	D_2 = G2.degree(node2)
	N_1 = G1.neighbors(node1)
	N_2 = G2.neighbors(node2)
	Sym_1 = G1.is_L_or_R(node1)
	Sym_2 = G2.is_L_or_R(node2)

	if Sym_1 == Sym_2:
		score = 0.5
	else:
		score = 0

	if not D_1 or not D_2:
		score = 0.5
		if Sym_1 == Sym_2:
			score += 0.25
		return score
	Aligned_neighbs = [ neighb for neighb in N_1 if neighb in aln and aln[neighb] in N_2 ]
	score += len( Aligned_neighbs ) / D_1
	return score

def extend_aln(G1, G2, aln, n1, n2, leftovers=False):
	""" attempt to add nodes to outer perimeter of subgraphs that 
	maintain isomorphism between these subgraphs and so building 
	the largest common subgraph """

	# for each of n1, align the node with the node in n2 with the highest similairty
	for node in n1:
		# build dictionary of similairites between node and all available 
		# nodes in the opposing graph
		if aln[node]:
			continue
		Similarities = { node2: Node_sim(G1, G2, node, node2, aln) for node2 in n2 }

		if Similarities:
			# Assign alignment of first best node based on similarity scores
			aln[node] = max(Similarities, key=Similarities.get)

			# remove aligned nodes from subset that are unaligned
			n2.remove(aln[node])
			
			if G1.is_L_or_R(node) and G2.is_L_or_R(aln[node]) and not leftovers:
				reflection = G1.reflect_n(aln[node])
				if reflection in n2:
					aln[G1.reflect_n(node)] = reflection
					n2.remove(reflection)
				
	return aln

def complete_ran_align(G1, G2, source=[["body", "body"]]):
	""" Randomly realign some nodes and check for improvement """
	# start alignment with the known node correspondance(s) from source
	aln = dict.fromkeys(G1.nodes())
	for i in source:
		aln[i[0]] = i[1]
	
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

	aln = extend_aln(G1, G2, aln, n1, n2, leftovers=True)

	# For complete data choose alignment with best symmetric substructure score
	return S3(G1, G2, aln), aln

def incomplete_ran_align(G1, G2, source=[["body", "body"]]):
	""" Randomly realign some nodes and check for improvement """
	# start alignment with the known node correspondance(s) from source
	aln = dict.fromkeys(G1.nodes())
	for i in source:
		aln[i[0]] = i[1]
	
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

	aln = extend_aln(G1, G2, aln, n1, n2, leftovers=True)

	# For complete data choose alignment with best symmetric substructure score
	return ics_score(G1, G2, aln), aln

def GrEvAl(Graph1, Graph2, repeat = 10, source=[["body", "body"]]):
	""" Grow an alignment from source, 
	the known node correspondance such as the body """

	G1 = Graph1.copy()
	G2 = Graph2.copy()
	best_score = 0
	best_aln = None
	if G1.number_of_nodes() > G2.number_of_nodes():
		G1, G2 = G2, G1
	

	if G1.graph['completeness'] + G2.graph['completeness'] == 2:
		# Make the graphs the same size so that a global alignment can be made
		if G1.number_of_nodes() > G2.number_of_nodes():
			G1, G2 = G2, G1
		
		
		for i in range(repeat):
			score, aln = complete_ran_align(G1.copy(), G2.copy(), source)
			if score > best_score:
				best_score = score
				best_aln = aln
		
		return best_score, best_aln
	
	else:
		
		for i in range(repeat):
			G1, G2 = G2, G1
			score, aln = incomplete_ran_align(G1.copy(), G2.copy(), source)
			if score > best_score:
				best_score = score
				best_aln = aln
		
		return best_score, best_aln



### Main ###
def main(argv):
	G1 = Plates.from_edgelist(argv[1])
	G2 = Plates.from_edgelist(argv[2])

	score, aln = GrEvAl(G1, G2, repeat=int(argv[3]))
	
	pprint.pprint(aln)
	print(score)
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)