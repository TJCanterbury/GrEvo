#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measurer distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'NewGrEvo.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from Glob import GrEvAl
from Morphlings import Plates
import numpy as np
import random

### Evo (hill climb) algorithm ###
def measurer(G1, G2):
	""" if isomorphic score is 0, else use S3 from a GrEvAl alignment """
	score, aln = GrEvAl(G1, G2)
	if nx.is_isomorphic(G1, G2):
		return 0, aln
	score = 1 - score
	
	return score, aln

class mutant():
	""" Morphling moves for a network of tectonic plate like characters, 
	with physical contact between characters giving an edge """
	def __init__(self, Score, Graph, move, aln):
		self.Score = Score
		self.Graph = Graph
		self.Aln = aln
		self.Move_dict = {
			0:"char_grows",
			1:"char_shrinks",
			2:"char_moves",
			3:"char_gain",
			4:"char_loss",
			5:"char_merge",
			6:"char_split",
			7:"char_expansion",
			8:"char_squeein",
			9:"char_squeeout"
		}
		self.Move = self.Move_dict[move]
	
	@classmethod
	def from_SnT(cls, G1, G2):
		Graph, Move = G1.mutator()
		while nx.is_isomorphic(Graph, G1):
			Graph, Move = G1.mutator()

		Score, Aln = measurer(Graph, G2)
		return cls(Score, Graph, Move, Aln)

	def __str__(self):
		self.Graph.__str__()
		return str(self.Aln) +"\n"+ str(self.Score)+"\n"+ \
			str(self.Move)+"\n"

def char_par(G1, G2, node1, node2):
	""" Biological distance between characters """		
	pars = 0
	if node1 and node2:
		for i in set(G1.nodes[node1]) | set(G2.nodes[node2]):
			pars += G1.nodes[node1][i] != G2.nodes[node2][i]

	return pars

def total_char_par(aln, G1, G2):
	""" Count number of character differences between aligned nodes """
	pars = 0
	
	if G1.graph["Align"]:
		for a in aln:
			pars += char_par(G1, G2, aln[a], a)
	else:
		for a in aln:
			pars += char_par(G2, G1, aln[a], a)

	return pars

def SA_GrEv(G1, G2, G1_name="a", G2_name="b", goal=100, temp = 0.1):
	""" finds the best next move for a given Morphling """
	Generation = 0
	best_score, old_aln = measurer(G1, G2)
	parsimony = 0
	G1_reset = G1.copy()
	reset = 2000
	
	while best_score != 0 and parsimony <= goal:
		# If dead-end reached reset
		if Generation > reset:
			print(G1_name + " X " + G2_name + ": " + "Failed, but will try again. Got stuck at: parsimony = " \
				 + str(parsimony) + ", Generation = " + str(Generation) + ", Score = " + str(best_score))
			reset *= 2
			Generation = 0
			G1 = G1_reset.copy()
			best_score, old_aln = measurer(G1, G2)
			parsimony = 0

		# Make a random move and measurer the effect
		M = mutant.from_SnT(G1=G1, G2=G2)
		Generation += 1

		if nx.is_isomorphic(M.Graph, G2) or M.Score == 0:
			best_M = M
			best_score = best_M.Score
			G1 = M.Graph
			parsimony += 1
			print(best_M)
			print(G1_name + " X " + G2_name + ": " + str(parsimony))
			break
		
		diff = M.Score - best_score
		t = temp / float(Generation + 1)
		metropolis = np.exp(-diff / t)
		
		# New state becoems current state if score improved or by chance depending on cost
		if M.Score < best_score and random.random() < metropolis and M.Score != best_score:
			best_M = M
			best_score = best_M.Score
			G1 = M.Graph
			print(best_M)
			parsimony += 1
	
	#parsimony += total_char_par(best_M.Aln, best_M.Graph, G2)

	return Generation, parsimony

### Business End ###
def main(argv):

	G1 = Plates.from_edgelist(argv[1])
	G2 = Plates.from_edgelist(argv[2])
	G1.graph['completeness'] = 1
	G2.graph['completeness'] = 1

	Generation, parsimony = SA_GrEv(G1, G2, temp=float(argv[3]))
	print(Generation)
	print(parsimony)
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)