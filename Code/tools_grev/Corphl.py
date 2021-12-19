#!/usr/bin/env python3

""" The Morphling and daughter classes for evolving from one species to another """

__appname__ = 'Morphlings.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
from graph_tool.all import *
import csv
import pandas as pd
from numpy import genfromtxt

## Classes ##
class Morph(Graph):
	""" Adds useful functions to the graph_tool graph class """
	def __init__(self, directed=False, prune=False, vorder=None):
		Graph.__init__(self, directed=False, prune=False, vorder=None)
		self.set_fast_edge_removal(fast=True)

	@classmethod
	def from_edgelist(cls, dir):
		""" derives Morph from edgelist format """
		g = cls()
		if 1 > 2:
			nodes = {}
			attrs = pd.read_csv(dir + "/C_Data.txt", sep=" ", index_col=0)
			name = g.new_vp("string")
			size = g.new_vp("short")
			tissue = g.new_vp("short")
			homolog = g.new_vp("short")
			for index, row in attrs.iterrows():
				Node = g.add_vertex()
				
		G_Data = genfromtxt(dir + "G_Data.txt", dtype=int)
		
		g.vp.id = g.add_edge_list(G_Data, hashed=True)
		return g 
	
	def draw(self, file="display.png"):
		""" Draw network, save to png"""
		graph_draw(self, output = file, vertex_text=self.vp.id, \
		 vertex_font_size=10, bg_color="white")



class biMorph(Morph):
	""" Adds useful functions to the graph_tool graph class """
	def __init__(self, directed=False, prune=False, vorder=None):
		Morph.__init__(self, directed=False, prune=False, vorder=None)

	def Set_Node_Symmetry(self, nodeInd):
		""" Assign symmetry VertexProperty values to node """
		node = self.vp.id[nodeInd]
		self.vp.Symmetry[nodeInd] = node % 10
	
	def Initialise_Symmetry(self):
		""" Set symmetry for all nodes """
		Symmetry = self.new_vp("int")
		self.vp.Symmetry=Symmetry
		for v in self.vertices():
			self.Set_Node_Symmetry(v)
	
	def DetectNeighbSymmetry(self, nodeInd):
		""" Detect Symmetry based on neighbour symmetry """
		sym = 0
		for w in nodeInd.all_neighbours():
			if self.vp.Symmetry[w] == 1:
				sym -= 1
			elif self.vp.Symmetry[w] == 2:
				sym += 1
		if sym > 0:
			return 2
		if sym < 0: 
			return 1
		return 0
	
	def ReflectNodeId(self, node):
		""" return id of nodes reflection """
		nId = self.vp.id[node]
		sym = nId % 10
		if sym > 0:
			if sym > 1:
				antId = nId - 1
			else:
				antId = nId + 1
		else:
			antId = nId
		return find_vertex(self, self.vp.id, antId)

	def ReflectNodes(self, nodes):
		""" return the reflected node indices """
		antinodes = []
		for v in nodes:
			antinodes.append(self.ReflectNodeId(v))
		return antinodes

	def NodeCluster(self, node):
		""" return node and it's neighbours """
		return list(self.get_all_neighbors(node)).append(node)

	def ExpectedClusReflection(self, node):
		""" Return hypothesized reflection of node and it's neighbours"""
		anticluster = self.ReflectNodes(NodeCluster(node))
		return anticluster

	def TestNodeTopoSymmetry(self, node):
		""" Test if reflection of node has same reflected
		neighbours """
		antinode = self.ReflectNodes(node)
		if set(self.NodeCluster(antinode)) == set(self.ExpectedClusReflection(node)):
			return True
		else:
			return False

	def FixNodeSymmetryValue(self, nodeInd):
		""" Set or correct Symmetry of a node """
		self.vp.Symmetry[nodeInd] = self.DetectNeighbSymmetry(nodeInd)
	
	def FixNodeSymTopo(self, node):
		""" Correct the topological symmetry for 
		symmetric property nodes """
		antinode = self.ReflectNodeId(node)
		TrueReflection = set(self.NodeCluster(antinode))
		ExpectedReflection = set(self.ExpectedClusReflection(node))
		if TrueReflection == ExpectedReflection:
			return
		else:
			for v in self.all_neighbors(node):
				self.remove_edge((node, v))
			for w in self.all_neigbors(antinode):
				self.add_edge(node, w, add_missing=True)

	def FixMorphSymmetry(self):
		""" fix symmetry value and topology for 
		all nodes, choosing left nodes as truth"""
		lefties = find_vertex(self, self.vp.Symmetry, 1)
		for v in lefties:
			self.FixNodeSymTopo(v)

## Functions ##
## Main ##
def main(argv):
	""" Test mutator """	
	G1 = biMorph.from_edgelist(argv[1])
	G1.Initialise_Symmetry()
	G1.FixMorphSymmetry()
	G1.draw()
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)