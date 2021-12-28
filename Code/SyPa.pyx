#!/usr/bin/env python3

""" Align 2 networks with the new alignment algortihm SYPA 
(Symmetrically aligned anatomical paths) and new alignment score CSS (Complete Symmetry Score) """

__appalnme__ = 'SyPa2.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from Morphlings import Placoderm
import random
import numpy as np

#### Functional Approach (currently 20 fold faster) ####
### Alignment Scores ###
def CSS(G1, G2, aln):
	""" Calculate symmetric substructure score but including all edges of target graph in denominator
	instead of aligned subgraph, with our goal being isomorphism instead of embedding of the source graph
	into the target. """
	edges_1 = set(G1.edges())
	edges_2 = set(G2.edges())
	c: cython.int = 0

	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]) or G2.has_edge(aln[e[1]], aln[e[0]]):
			c += 1

	score = c / ( len(edges_1) - c + len(edges_2))
	
	return score

def S3(G1, G2, aln):
	""" Calculate symmetric substructure score described by Saraph and Milenković
	(Saraph V. Milenković T. (2014) MAGNA: maximizing accuracy in global network alignment. Bioinformatics, 30, 2931–2940.) """
	edges_1 = set(G1.edges())
	
	c: cython.int = 0
	subnodes = []

	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]) or G2.has_edge(aln[e[1]], aln[e[0]]) \
			and G1.is_L_or_R(e[0]) == G1.is_L_or_R(aln(e[0])) \
			and G1.is_L_or_R(e[1]) == G1.is_L_or_R(aln(e[1])):
			c += 1
			subnodes.append(aln[e[0]])
			subnodes.append(aln[e[1]])

	subG1 = G2.subgraph(subnodes)
	sublen: cython.int = subG1.size()
	score = c /( sublen + len(edges_1) - c)
	return score

### Node similarity scores ###
def Char_sim(G1, G2, node1, node2, D1, D2):
	""" Biological distance between characters """
	if not max(len(G1.nodes[node1]),len(G2.nodes[node2])) :
		return 1
		
	pars: cython.int = 0
	for i in set(G1.nodes[node1]).union(set(G2.nodes[node2])):
		if i in G1.nodes[node1] and i in G2.nodes[node2]:
			pars += G1.nodes[node1][i] != G2.nodes[node2][i]

	pars += D1 != D2

	score = pars / (max(len(G1.nodes[node1]),len(G2.nodes[node2])) + 1)
	score = 1 - score
	return score 

def Node_sim(G1, G2, node1, node2, aln):
	""" Determine node similairty based on symmetry,
	and alignment """
	score: cython.int

	try:
		Sym_1: cython.int = node1 % 10
		Sym_2: cython.int = node2 % 10
	except:
		print(node1)
		print(G1.nodes.data())
		print(node2)
		print(G2.nodes.data())
		exit

	if Sym_1 == Sym_2:
		score = 0.5
	else:
		return 0
	
	N_1: cython.list = G1._adj[node1]
	N_2: cython.list = G2._adj[node2]
	D_1: cython.int = len(N_1)

	Aligned_neighbs: cython.int = 0
	for neighb in N_1:
		if neighb in aln and aln[neighb] in N_2:
			Aligned_neighbs += 1 
	score += Aligned_neighbs / D_1
	# uncomment below for use of character info:
	#D_2: cython.int = len(N_2)
	#bio_score: cython.int = Char_sim(G1, G2, node1, node2, D_1, D_2)
	#score = (score + bio_score) / 2
	return score

### Align Homologues ###
def User_aln(G1, G2, aln):
	""" Align nodes with user defined homology 
	under the homolog column of the Character data, 
	provided and encoded as node attributes """
	ua1 = nx.get_node_attributes(G1, "homolog")
	ua1 = {x:y for x,y in ua1.items() if y!=0}
	ua2 = nx.get_node_attributes(G2, "homolog")
	ua2 = {x:y for x,y in ua2.items() if y!=0}

	for node in ua1:
		n2 = [n for n in ua2 if ua2[n] == ua1[node]]
		aln = Sym_align_nodes(G1, G2, aln, [node], n2)
	return aln

### Align Given nodes ###
def Sym_align_nodes(G1, G2, aln, n1, n2):
	""" Align nodes that have the maximum node similarity, could be improved with the Hungarian algorithm
	if not for coneighbours being a factor in node similarity. Change the node similarity function to use 
	biological and graphlet degree similarity and use of hungarian algorithm would be optimal. """
	
	n1 = [n for n in n1 if not aln[n]]
	n2 = [n for n in n2 if n not in aln.values()]
	random.shuffle(n1)

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
			
			reflection1 = G1.Brother_V(node)
			reflection2 = G2.Brother_V(aln[node])
			if reflection1 in n1 and reflection2 in n2 and node%10 and aln[node]%10:
				aln[reflection1] = reflection2
				n2.remove(reflection2)
	
	return aln

### Unalign central nodes from symmetric nodes ###
def Correct_sym(aln):
	""" Unalign nodes that are of a different symmetry attribute """
	for k in aln:
		try:
			if k%10 != aln[k]%10:
				aln[k] = None
		except:
			pass
	return aln

def Grow_aln(aln, G1, G2):
	""" Greedily extend alignment, to 1 or 2 radii out at a time """
	
	# Pick random source:
	source = []
	homologs = list({ele for ele in aln if aln[ele]})
	if len(homologs) > 1:
		source.append( random.choice(homologs) )
	else:
		source.append(homologs[0])
	source.append(aln[source[0]])
	
	# make independent alignment from each node in source
	n1 = nx.single_source_shortest_path_length(G1, source[0])
	n2 = nx.single_source_shortest_path_length(G2, source[1])
	max_radius = min([max(n1.values()), max(n2.values())])
	radius = 1

	while radius <= max_radius:
		A_neighbs = [k for k in n1 if n1[k] == radius] 
		B_neighbs = [k for k in n2 if n2[k] == radius] 
		aln = Sym_align_nodes(G1, G2, aln, A_neighbs, B_neighbs)
		radius += np.random.choice([1,2])
	
	return aln

### User align, greedily extend alignment and correct symmetry
def Align(G1, G2, aln, Complete = True):
	"""  	
		Greedily extend alignment, fill in the 
		gaps and correct for symmetry. 
	"""

	aln = Grow_aln(aln, G1, G2)

	#Align any remaining nodes
	aln = Sym_align_nodes(G1, G2, aln, G1.nodes(), G2.nodes())
	aln = Correct_sym(aln)

	if Complete:
		# For complete data choose alignment with best symmetric substructure score
		return CSS(G1, G2, aln), aln
		
	else:
		# For incomplete data choose alignment with best symmetric substructure score
		# so as to obtimise for embedding instead of matching the source to the target
		return S3(G1, G2, aln), aln

def Improve(aln, G1, G2, score1, repeats = 1, Complete=True):
	""" Identify nodes that aligned poorly and test different alignments"""
	
	edges_1 = set(G1.edges())
	bad_nodes = []

	for e in edges_1:
		if not G2.has_edge(aln[e[0]], aln[e[1]]) or not G2.has_edge(aln[e[1]], aln[e[0]]):
			bad_nodes.append(e[0])
			bad_nodes.append(e[1])
	
	good_nodes = set(G1.nodes()) - set(bad_nodes)
	aln0 = dict.fromkeys(G1.nodes())
	aln0 = User_aln(G1, G2, aln0)

	# Delete alignments that are not good:
	for i in good_nodes:
		aln0[i]= aln[i]

	#Try to improve it
	for i in range(repeats):
		aln2 = Sym_align_nodes(G1, G2, aln0.copy(), G1.nodes(), G2.nodes())
		aln2 = Correct_sym(aln2)

	if Complete:
		score2 = CSS(G1, G2, aln2)
	else:
		score2 = S3(G1, G2, aln2)

	if score2 > score1:
		return score2, aln2
	else:
		return score1, aln

### Tool to check that a global 1 to 1 alignment was generated
def Check_global(aln):
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

### Repeat alignment
def SYPA(G1, G2, repeat = 25):
	""" Grow an alignment from source, 
	the known node correspondance such as the body """
	if G1.graph['completeness'] + G2.graph['completeness'] == 2:
		Complete = True
	else:
		Complete = False
	
	best_score = 0
	best_aln = None
	
	aln0 = dict.fromkeys(G2.nodes())
	
	#set user defined alignments
	aln0 = User_aln(G2, G1, aln0)
	aln0 = Correct_sym(aln0)

	for i in range(repeat):
		score, aln = Align(G2, G1, aln=aln0.copy(), Complete=Complete)
		#score, aln = Improve(aln, G2, G1, score1=score, Complete=Complete)

		if score > best_score:
			best_score = score
			best_aln = aln
			G1.graph["Align"] = False
			G2.graph["Align"] = True

	return best_score, best_aln

### Main ###
def main(argv):	
	G1 = Placoderm.From_Dir(argv[1])	
	G2 = Placoderm.From_Dir(argv[2])	
	
	score, aln = SYPA(G1, G2, repeat=int(argv[3]))		
	
	print(aln)	
	print(score)	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)