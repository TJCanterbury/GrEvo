#!/usr/bin/env python3

""" Align 2 networks with the new alignment algortihm SYPA 
(Symmetrically aligned anatomical paths) and new alignment score CSS (Complete Symmetry Score) """

__appalnme__ = 'SYPA.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from Morphlings import Plates
import random
import pprint

### Align ###
def CSS(G1, G2, aln):
	""" Calculate symmetric substructure score but including all edges of target graph in denominator
	instead of aligned subgraph, with our goal being isomorphism instead of embedding of the source graph
	into the target. """
	edges_1 = set(G1.edges())
	edges_2 = set(G2.edges())
	c = 0

	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]) or G2.has_edge(aln[e[1]], aln[e[0]]):
			c += 1

	score = c / ( len(edges_1) - c + len(edges_2))
	
	return score

def S3(G1, G2, aln):
	""" Calculate symmetric substructure score described by Saraph and Milenković
	(Saraph V. Milenković T. (2014) MAGNA: maximizing accuracy in global network alignment. Bioinformatics, 30, 2931–2940.) """
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

def user_aln(G1, G2, aln):
	""" Align nodes with user defined homology 
	under the homolog column of the Character data, 
	provided and encoded as node attributes """

	ua1 = nx.get_node_attributes(G1, "homolog")
	ua1 = {x:y for x,y in ua1.items() if y!=0}
	ua2 = nx.get_node_attributes(G2, "homolog")
	ua2 = {x:y for x,y in ua2.items() if y!=0}

	for node in ua1:
		n2 = [n for n in ua2 if ua2[n] == ua1[node]]
		aln = extend_aln(G1, G2, aln, [node], n2)
	return aln

def char_sim(G1, G2, node1, node2):
	""" Biological distance between characters """
	if not max(len(G1.nodes[node1]),len(G2.nodes[node2])) :
		return 1
		
	pars = 0
	for i in set(G1.nodes[node1]).union(set(G2.nodes[node2])):
		if i in G1.nodes[node1] and i in G2.nodes[node2]:
			pars += G1.nodes[node1][i] != G2.nodes[node2][i]

	pars += G1.degree(node1) != G2.degree(node2)

	score = pars / (max(len(G1.nodes[node1]),len(G2.nodes[node2])) + 1)
	score = 1 - score
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
	bio_score = char_sim(G1, G2, node1, node2)
	score = (score + bio_score) / 2
	return score

def extend_aln(G1, G2, aln, n1, n2):
	""" Align nodes that have the maximum node similarity, could be improved with the Hungarian algorithm
	if not for coneighbours being a factor in node similarity. Change the node similarity function to use 
	biological and graphlet degree similarity and use of hungarian algorithm would be optimal. """
	n1 = [n for n in n1 if not aln[n]]
	n2 = [n for n in n2 if n not in aln.values()]
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
			
			reflection1 = G1.reflect_n(node)
			reflection2 = G1.reflect_n(aln[node])
			if reflection1 in n1 and reflection2 in n2 and G1.is_L_or_R(node) and G2.is_L_or_R(aln[node]):
				aln[reflection1] = reflection2
				n2.remove(reflection2)
				
	return aln

def correct_sym(aln, G1, G2):
	""" Unalign nodes that are of a different symmetry attribute """
	for k in aln:
		try:
			if G1.is_L_or_R(k) != G1.is_L_or_R(aln[k]):
				aln[k] = None
		except:
			pass
	return aln

def ran_align(G1, G2, source=[["body", "body"]]):
	""" Randomly realign some nodes and check for improvement """
	# start alignment with the known node correspondance(s) from source
	aln = dict.fromkeys(G1.nodes())
	for i in source:
		aln[i[0]] = i[1]
	
	#set user defined alignments
	aln = user_aln(G1, G2, aln)

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
	
	#Align any remaining nodes
	aln = extend_aln(G1, G2, aln, G1.nodes(), G2.nodes())

	aln = correct_sym(aln, G1, G2)

	if G1.graph['completeness'] + G2.graph['completeness'] == 2:
		# For complete data choose alignment with best symmetric substructure score
		return CSS(G1, G2, aln), aln
	else:
		# For incomplete data choose alignment with best symmetric substructure score
		# so as to obtimise for embedding instead of matching the source to the target
		return S3(G1, G2, aln), aln

def check_global(aln):
	""" Check that an alignment gives a one to one mapping (very important it does) """
	rev_multidict = {}
	for key, value in aln.items():
		rev_multidict.setdefault(value, set()).add(key)
	keys = [key for key, values in rev_multidict.items() if len(values) > 1]
	if keys:
		print("Not a one to one mapping, the following nodes aligned more than once: ")
		print(keys)
		return 1
	
	return 0

def SYPA(G1, G2, repeat = 20, source=[["body", "body"]]):
	""" Grow an alignment from source, 
	the known node correspondance such as the body """

	G1.add_symmetry()
	G2.add_symmetry()
	best_score = 0
	best_aln = None

	for i in range(repeat):
		score, aln = ran_align(G2, G1, source)
		if score > best_score:
			best_score = score
			best_aln = aln
			G1.graph["Align"] = False
			G2.graph["Align"] = True
	
	return best_score, best_aln

### Main ###
def main(argv):
	G1 = Plates.from_edgelist(argv[1])
	G2 = Plates.from_edgelist(argv[2])
	G1.graph['completeness'] = int(argv[4])
	G2.graph['completeness'] = int(argv[5])

	score, aln = SYPA(G1, G2, repeat=int(argv[3]))
	
	pprint.pprint(aln)
	print(score)
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)