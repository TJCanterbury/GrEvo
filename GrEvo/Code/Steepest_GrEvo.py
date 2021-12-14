#!/usr/bin/env python3

""" Steepest ascent hill climb from one graph to another with bio/Morphology data """

__appname__ = 'Steepest_GrEvo.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
from Morphlings import Plates
import networkx as nx
from SYPA2 import SYPA
from Morphlings import Plates


### Evo (hill climb) algorithm ###
def measurer(G1, G2):
	""" if isomorphic score is 0, else use S3 from a SYPA alignment """
	score, aln = SYPA(G1, G2)
	if nx.is_isomorphic(G1, G2):
		0, aln
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

	def plot(self, Pars):
		file = str(Pars)+".png"
		self.Graph.__str__(file)
		print(str(self.Aln) +"\n"+ str(self.Score)+"\n"+ \
			str(self.Move)+"\n")
		return 0

def tissue_par(G1, G2, node1, node2):
	""" Tissue difference between characters """		
	pars = 0
	i = "tissue_type"
	
	if not node1 or not node2:
		return 1

	if i in G1.nodes[node1] and i in G2.nodes[node2]:
		pars += G1.nodes[node1][i] != G2.nodes[node2][i]

	else:
		if i in G1.nodes[node1]:
			G2.nodes[node2][i] = G1.nodes[node1][i]
		if i in G2.nodes[node2]:
			G1.nodes[node1][i] = G2.nodes[node2][i]

	return pars

def total_char_par(aln, G1, G2):
	""" Count number of character differences between aligned nodes """
	pars = 0

	if G1.graph["Align"]:
		for a in aln:
			pars += tissue_par(G1, G2, a, aln[a])
	else:
		for a in aln:
			pars += tissue_par(G2, G1, a, aln[a])

	return pars

def Steep_GrEv(G1, G2, G1_name="a", G2_name="b", goal=100, Breadth=10000, printer = False): 
	""" finds the best next move for a given Morphling """
	Generation = 0
	best_score, old_aln = measurer(G1, G2)
	parsimony = 0
	reset = 20000
	G1_reset = G1.copy()
	improved=False

	if printer:
		print(old_aln)
		print(best_score)
		print("\n")

	while best_score != 0 and parsimony <= goal:
		# Make a random move and measurer the effect
		M = mutant.from_SnT(G1=G1, G2=G2)
		Generation += 1

		if M.Score == 0:
			best_M = M
			best_score = best_M.Score
			G1 = M.Graph
			parsimony += 1
			if printer:
				best_M.plot(Pars=parsimony)
			break
		
		# New state becomes current state if score improved
		if M.Score < best_score and M.Score != best_score:
			best_M = M
			best_score = best_M.Score
			improved = True
		
		if Generation >= Breadth and improved:
			G1 = best_M.Graph
			parsimony += 1
			Generation = 0
			improved=False
			if printer:
				best_M.plot(Pars=parsimony)
		
		if Generation > reset:
			print(G1_name + " X " + G2_name + ": " + "Failed, but will try again. Got stuck at: parsimony = " \
				 + str(parsimony) + ", Generation = " + str(Generation) + ", Score = " + str(best_score))
			Generation = 0
			parsimony = 0
			reset *= 2
			G1 = G1_reset.copy()
			best_score, old_aln = measurer(G1, G2)

	parsimony += total_char_par(best_M.Aln, best_M.Graph, G2)

	return Generation, parsimony

### Business End ###
### Business End ###
def main(argv):
	G1 = Plates.from_edgelist(argv[1]+"G_Data.txt")
	G2 = Plates.from_edgelist(argv[2]+"G_Data.txt")
	G1.graph['completeness'] = 1
	G2.graph['completeness'] = 1
	G1.attr_from_csv(argv[1]+"C_Data.txt")
	G2.attr_from_csv(argv[2]+"C_Data.txt")

	Generation, parsimony = Steep_GrEv(G1, G2, goal=float(argv[3]), Breadth=float(argv[4]), printer=True)
	print(Generation)
	print(parsimony)
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)