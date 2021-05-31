#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measurer distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'GrEvo_Synth.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.4'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import networkx as nx
from nxcode import readx
from nxcode import drawx
import gressure as gr

## Classes ##
class Morphling:

	def __init__(self, Graph):
		self.G = Graph
		self.degree = self.G.degree()
		self.nodes = self.G.nodes()
		self.number_of_nodes = self.G.number_of_nodes()

	def Add_edge(self, Node1, Node2):
		reflection = gr.Reflect([Node1, Node2])
		self.G.add_edge(Node1, Node2)
		self.G.add_edge(reflection[0], reflection[1])
	
	def Remove_node(self, Node):
		reflection = gr.Reflect([Node])
		self.G.remove_nodes_from((Node, reflection[0]))

	def Mean_degree(self):
		my_degrees = self.G.degree()
		degree_values = [v for k, v in my_degrees]
		sum_G = sum(degree_values)

		mean = (sum_G / self.G.number_of_nodes()) 
		self.mean_degree = mean
		return mean

## Functions ##
### General Graph Functions ###
def num_edges(G):
	""" Returns from a normal distribution an integer for the
	number of edges to add """
	
	randomInts = np.random.normal(loc=G.mean_degree, size=1).astype(int)
	
	while randomInts < 0:
		randomInts = np.random.normal(loc=G.mean_degree, size=1).astype(int)
	
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

### Business End ###
def main(argv):
	G1 = readx(argv[1])
	G1 = Morphling(G1)
	G1.Add_edge("L9", "L90")
	print(G1.nodes)
	print(G1.Mean_degree())
	print(G1.degree)
	
	drawx(G1.G, fname="display_test.png")
	
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)