#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measurer distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'GrEvo.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.5'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from Glob import GrEvAl
from Morphlings import Plates

### Evo (hill climb) algorithm ###
def measurer(G1, G2, ret_aln = False):
	""" if isomorphic score is 0, else use S3 from a GrEvAl alignment """
	score, aln = GrEvAl(G1.copy(), G2.copy())

	score = 1 - score
	if ret_aln:
		return score, aln
	return score

def recorder(move, score, print_ = True):
	""" Store moves that improve portrait divergence score for future optimisation """

	Line = "Move: " + str(move) + " Score: " + str(score)+"\n"

	if print_ == True:
		print(Line)
	return Line

def searcher(G1, G2, size, attempts, dead_ends, old_score):
	""" finds the best next move for a given Morphling """
	breadth = 0
	best_G = 0
	best_move = None
	best_score = float(1)
	stuck = 0
	trying = 0

	while breadth < size or best_score >= old_score:
		# Make a random move and measurer the effect
		morph, move = G1.mutator()
		
		score, aln = measurer(G1=morph, G2=G2, ret_aln=True)
		breadth += 1

		if score == 0 or nx.is_isomorphic(morph, G2):
			return 0, move, morph, aln

		# Escape dead end and record its graph
		if stuck >= attempts:    
			return old_score, None, G1, None
		
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
		if score < best_score and not move == None:
			best_score = score
			best_G = morph
			best_move = move
			best_aln = aln

		# Record how stuck we are on this generation
		if breadth > size:
			stuck += 1
			breadth = 0
			print(stuck)
	print(best_aln)
	return best_score, best_move, best_G, best_aln

def climber(G1, G2, sample_size=1000, attempts=3, goal=100):
	""" Apply mutator, if score improved record and recurse """
	Generation = 0
	moves1 = []
	g1s = [G1]
	old_score = measurer(G1, G2)
	move = "Start Distance"
	moves1.append("G1 " + recorder(move, old_score))
	dead_ends = []
	size = sample_size
	best_score = old_score
	
	while best_score != 0:
		if Generation + 1 == goal:
			dead_ends.append(best_G) 
			Generation = 0
			g1s = [G1]
			move = "Start Distance"
			moves1.append("G1 " + recorder(move, old_score))
			size = sample_size
			best_score = old_score

		best_score, best_move, best_G, best_aln = searcher(G1=g1s[Generation], 
			G2=G2, size=size, attempts=attempts, dead_ends=dead_ends,
			 old_score = best_score)
		
		if best_aln:
			Generation += 1
			size = sample_size * (Generation + 1)
			best_G.__str__()
			g1s.append(best_G.copy())
			moves1.append(recorder(best_move, best_score))
		
		elif len(dead_ends) >= 3:
			return Generation * 10

		else:
			dead_ends.append(best_G) 
			Generation = 0
			g1s = [G1]
			move = "Start Distance"
			moves1.append("G1 " + recorder(move, old_score))
			size = sample_size
			best_score = old_score
	
	return Generation

### Business End ###
def main(argv):

	if argv[3] == "manual":
		G1 = Plates.from_edgelist(argv[1])
		G2 = Plates.from_edgelist(argv[2])
		G1.graph['completeness'] = 1
		G2.graph['completeness'] = 1
		
		
		if len(argv) == 6:
			move=int(argv[5])
		else:
			move = None
		morph, move = G1.mutator(Node=argv[4], move=move)
		score = measurer(morph, G2)
		recorder(move, score)

		morph.__str__("display_test.png")
	
	else:
		G1 = Plates.from_edgelist(argv[1])
		G2 = Plates.from_edgelist(argv[2])
		G1.graph['completeness'] = 1
		G2.graph['completeness'] = 1

		if G1.graph['completeness'] < G2.graph['completeness']:
			print("Evolving from " + argv[2] + " due to it's greater completeness")
			G1, G2 = G2, G1

		parsimony = climber(G1, G2, sample_size=int(argv[3]))

		print(parsimony)
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)