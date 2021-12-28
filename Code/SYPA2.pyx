#!/usr/bin/env python3

""" Align 2 networks with the new alignment algortihm SYPA 
(Symmetrically aligned anatomical paths) and new alignment score CSS (Complete Symmetry Score) """

__appalnme__ = 'SyPa3.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from Morphlings import Placoderm
import random
import numpy as np

#### Oop approach ####
class Alignment():
	def __init__(self, G1, G2):
		self.aln = dict.fromkeys(G1.nodes())
		self.source = G1
		self.target = G2
		self.best_score = 0
		self.score = 0
		self.s_num_e = G1.number_of_edges()
		self.t_num_e = G2.number_of_edges()
		self.s_num_v = G1.number_of_nodes()
		self.t_num_v = G2.number_of_nodes()

	@property
	def score(self):
		return self.score
	
	@score.setter
	def score(self, value): # called whenever self.score is changed (in CSS and S3)
							# to set the new best score and alignment
		if value > self.best_score:
			self.best_score = value
			self.best_aln = self.aln.copy()
	
	def align_node_pair(self, u, v):
		self.aln[u] = v
		self.aln[self.target.Brother_V(u)] = self.target.Brother_V(v)
	
	def count_edge_match(self):
		c = 0
		for e in self.source.edges():
			if self.target.has_edge(self.aln[e[0]], self.aln[e[1]]) or self.target.has_edge(self.aln[e[1]], self.aln[e[0]]):
				c += 1
		return c
	
	def CSS(self):
		""" Calculate symmetric substructure score but including all edges of target graph in denominator
		instead of aligned subgraph, with our goal being isomorphism instead of embedding of the source graph
		into the target. """
		c = self.count_edge_match()

		self.score = c / ( self.s_num_e - c + self.t_num_e)
	
	def S3(self):
		""" Calculate symmetric substructure score described by Saraph and Milenković
		(Saraph V. Milenković T. (2014) MAGNA: maximizing accuracy in global network alignment. Bioinformatics, 30, 2931–2940.) """
		edges_1 = set(self.source.edges())
		G1 = self.source
		G2 = self.target

		c = 0
		subnodes = []

		for e in edges_1:
			if G2.has_edge(self.aln[e[0]], self.aln[e[1]]) or G2.has_edge(self.aln[e[1]], self.aln[e[0]]) \
				and e[0] % 10 == self.aln(e[0]) % 10 and e[1] % 10 == self.aln(e[1]) % 10 :
				c += 1
				subnodes.append(self.aln[e[0]])
				subnodes.append(self.aln[e[1]])

		subG1 = G2.subgraph(subnodes)
		sublen = subG1.size()
		self.score = c /( sublen + len(edges_1) - c)

	def Node_sim(self, u, v):
		""" Determine node similairty based on symmetry,
		and alignment """

		score: cython.int
		Sym_1: cython.int = u%10
		Sym_2: cython.int = v%10

		if Sym_1 == Sym_2:
			score = 0.5
		else:
			return 0
		
		N_1: cython.list = self.source._adj[u]
		N_2: cython.list = self.target._adj[v]
		D_1: cython.int = len(N_1)
		
		Aligned_neighbs: cython.int = 0
		for neighb in N_1:
			if neighb in self.aln and self.aln[neighb] in N_2:
				Aligned_neighbs += 1 
		score += Aligned_neighbs / D_1
		return score

	### Align Homologues ###
	def User_aln(self):
		""" Align nodes with user defined homology 
		under the homolog column of the Character data, 
		provided and encoded as node attributes """
		G1 = self.source
		G2 = self.target

		ua1 = nx.get_node_attributes(G1, "homolog")
		ua1 = {x:y for x,y in ua1.items() if y!=0}
		ua2 = nx.get_node_attributes(G2, "homolog")
		ua2 = {x:y for x,y in ua2.items() if y!=0}

		for node in ua1:
			n2 = [n for n in ua2 if ua2[n] == ua1[node]]
			self.Sym_align_nodes([node], n2)

	### Align Given nodes ###
	def Sym_align_nodes(self, n1, n2):
		""" Align nodes that have the maximum node similarity, could be improved with the Hungarian algorithm
		if not for coneighbours being a factor in node similarity. Change the node similarity function to use 
		biological and graphlet degree similarity and use of hungarian algorithm would be optimal. """

		G1 = self.source
		G2 = self.target
		n1 = [n for n in n1 if not self.aln[n]]
		n2 = [n for n in n2 if n not in self.aln.values()]
		random.shuffle(n1)

		# for each of n1, align the node with the node in n2 with the highest similairty
		for node in n1:
			# build dictionary of similairites between node and all available 
			# nodes in the opposing graph
			if self.aln[node]:
				continue
			Similarities = { node2: self.Node_sim(node, node2) for node2 in n2 }

			if Similarities:
				# Assign alignment of first best node based on similarity scores
				self.aln[node] = max(Similarities, key=Similarities.get)

				# remove aligned nodes from subset that are unaligned
				n2.remove(self.aln[node])

				reflection1 = G1.Brother_V(node)
				reflection2 = G2.Brother_V(self.aln[node])
				if reflection1 in n1 and reflection2 in n2 and G1.nodes[node]['Sym'] \
					and G2.nodes[self.aln[node]]['Sym']:
					self.aln[reflection1] = reflection2
					n2.remove(reflection2)

	### Unalign central nodes from symmetric nodes ###
	def Correct_sym(self):
		""" Unalign nodes that are of a different symmetry attribute """
		
		G1 = self.source

		for k in self.aln:
			try:
				if G1.nodes[k]['Sym'] != G1.nodes[self.aln[k]]['Sym']:
					self.aln[k] = None
			except:
				pass

	def Grow_aln(self):
		""" Greedily extend alignment, to 1 or 2 radii out at a time """

		G1 = self.source
		G2 = self.target
		# Pick random source:
		source = []
		homologs = list({ele for ele in self.aln if self.aln[ele]})
		if len(homologs) > 1:
			source.append( random.choice(homologs) )
		else:
			source.append(homologs[0])
		source.append(self.aln[source[0]])

		# make independent alignment from each node in source
		n1 = nx.single_source_shortest_path_length(G1, source[0])
		n2 = nx.single_source_shortest_path_length(G2, source[1])
		max_radius = min([max(n1.values()), max(n2.values())])
		radius = 1

		while radius <= max_radius:
			A_neighbs = [k for k in n1 if n1[k] == radius] 
			B_neighbs = [k for k in n2 if n2[k] == radius] 
			self.Sym_align_nodes(A_neighbs, B_neighbs)
			radius += np.random.choice([1,2])

	### User align, greedily extend alignment and correct symmetry
	def Align(self, Complete = True):
		"""  	
			Greedily extend alignment, fill in the 
			gaps and correct for symmetry. 
		"""

		G1 = self.source
		G2 = self.target

		self.Grow_aln()

		#Align any remaining nodes
		self.Sym_align_nodes(G1.nodes(), G2.nodes())
		self.Correct_sym()

		if Complete:
			# For complete data choose alignment with best symmetric substructure score
			self.CSS()

		else:
			# For incomplete data choose alignment with best symmetric substructure score
			# so as to obtimise for embedding instead of matching the source to the target
			self.S3()

	def Improve(self, repeats = 1, Complete=True):
		""" Identify nodes that aligned poorly and test different alignments"""

		G1 = self.source
		G2 = self.target

		edges_1 = set(G1.edges())
		bad_nodes = []

		for e in edges_1:
			if not G2.has_edge(self.best_aln[e[0]], self.best_aln[e[1]]) or not G2.has_edge(self.best_aln[e[1]], self.best_aln[e[0]]):
				bad_nodes.append(e[0])
				bad_nodes.append(e[1])

		good_nodes = set(G1.nodes()) - set(bad_nodes)
		self.aln = dict.fromkeys(G1.nodes())
		self.User_aln()

		# Delete alignments that are not good:
		for i in good_nodes:
			self.aln[i]= self.best_aln[i]

		#Try to improve it
		for i in range(repeats):
			self.Sym_align_nodes(G1.nodes(), G2.nodes())
			self.Correct_sym()

		if Complete:
			self.CSS()
		else:
			self.S3()

	### Repeat alignment
	def SYPA_aln(self, repeat = 25):
		""" Grow an alignment from source, 
		the known node correspondance such as the body """
		if self.source.graph['completeness'] + \
			self.target.graph['completeness'] == 2:
			Complete = True
		else:
			Complete = False

		#set user defined alignments
		self.User_aln()
		self.Correct_sym()

		for i in range(repeat):
			self.Align(Complete=Complete)
			self.Improve(Complete=Complete)

		return self.best_score, self.best_aln

def SYPA(G1, G2, repeat = 25):
	G1.graph["Align"] = False
	G2.graph["Align"] = True
	A = Alignment(G1=G2, G2=G1)
	return A.SYPA_aln(repeat)

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