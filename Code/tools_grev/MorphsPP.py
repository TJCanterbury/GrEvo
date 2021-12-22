#!/usr/bin/env python3

""" The Morphling and daughter classes for evolving from one species to another """

__appname__ = 'Morphlings.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import graph_tool.all as gt
import pandas as pd
from numpy import genfromtxt
import cProfile
import numpy as np
import itertools as it

## Classes ##
class Morph(gt.Graph):
	""" Adds useful functions to the graph_tool graph class """
	def __init__(self, directed=False, prune=False, vorder=None):
		gt.Graph.__init__(self, directed=False, prune=False, vorder=None)
		self.set_fast_edge_removal(fast=True)
		self.graph={}

	@classmethod
	def From_Dir(cls, dir, Completeness = 1):
		G = Placoderm.from_edgelist(dir)
		G.attr_from_csv(dir)
		G.graph['completeness'] = Completeness
		G.graph['dir']=dir
		G.Initialise_Symmetry()
		return G

	@classmethod
	def from_edgelist(cls, dir):
		""" derives Morph from edgelist format """
		g = cls()
				
		G_Data = genfromtxt(dir + "G_Data.txt", dtype=int)
		
		g.vp.id = g.add_edge_list(G_Data, hashed=True)
		return g 
	
	def attr_from_csv(self, dir):
		""" import character data from csv, add to properties 
		and save property titles """
		attrs = pd.read_csv(dir + "/C_Data.txt", sep=" ", index_col=0)
		self.attr = attrs.to_dict(orient='index')
		return 0
	
	def draw(self, file="display.png"):
		""" Draw network, save to png"""
		gt.graph_draw(self, output = file, vertex_text=self.vp.id, \
		 vertex_font_size=10, bg_color="white")
	
	def Neighbours(self, node):
		return self.get_all_neighbors(node)
	
	def homologs(self):
		""" return number of unique homologies """
		ua = self.vp.homolog
		return set(ua.values())
	
	def degree(self, v):
		degree = self.get_total_degrees([v])
		return degree
	
	def degrees(self):
		return self.get_total_degrees(self.get_vertices())
	
	def ran_degree(self):
		""" Returns from a normal distribution an integer for the
		number of edges to add """
		randomInts = np.random.normal(loc=0, scale=3, size=1).astype(int)
		max_degree = np.max(self.degrees())
		while randomInts < 0 or randomInts > max_degree:
			randomInts = np.random.normal(loc=0, size=1).astype(int)
		
		return int(randomInts)

	def NodeCluster(self, node):
		""" return node and it's neighbours """
		Cluster =  list(self.get_all_neighbors(node))
		Cluster.append(node)

		return Cluster

	def Has_Edge(self, u, v):
		""" return 1 if edge exists in graph between u and v """
		return v in self.Neighbours(u)
	
	def coneighbours(self, nodes):
		""" returns set of common neighbours between 2 nodes """
		cneighbs = set(self.Neighbours(nodes[0]))
		for v in nodes:
			cneighbs = cneighbs & set(self.Neighbours(v))
		return cneighbs

class biMorph(Morph):
	""" Adds useful bi lateral symmetry functions to the graph_tool graph class """
	def __init__(self, directed=False, prune=False, vorder=None):
		Morph.__init__(self, directed=False, prune=False, vorder=None)

	def Set_Node_Symmetry(self, nodeInd):
		""" Assign symmetry VertexProperty values to node """
		node = self.vp.id[nodeInd]
		self.vp.Symmetry[nodeInd] = node % 10
		self.vp.pairing[nodeInd] = node - self.vp.Symmetry[nodeInd]
		self.vp.brother[nodeInd] = self.find_brother(nodeInd)
	
	def Initialise_Symmetry(self):
		""" Set symmetry for all nodes """
		Symmetry = self.new_vp("int")
		self.vp.Symmetry=Symmetry
		pairing = self.new_vp("int")
		self.vp.pairing=pairing
		brother = self.new_vp("int")
		self.vp.brother=brother
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
	
	def Brother_Node(self, node):
		""" return index of nodes reflection """
		return self.vp.brother[node]

	def find_brother(self, node):
		""" return index of nodes reflection """
		nId = self.vp.id[node]
		sym = nId % 10
		if sym > 0:
			if sym > 1:
				antId = nId - 1
			else:
				antId = nId + 1
			return gt.find_vertex(self, self.vp.id, antId)[0]
		else:
			return node
		
	def ReflectNodes(self, nodes):
		""" return the reflected node indices """
		antinodes = []
		for v in nodes:
			antinodes.append(self.Brother_Node(v))
		return antinodes

	def ExpectedClusReflection(self, node):
		""" Return hypothesized reflection of node and it's neighbours"""
		anticluster = self.ReflectNodes(self.NodeCluster(node))
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

	def FixEdgeSymTopo(self, node):
		""" Correct the edge topology symmetry for 
		symmetric property nodes """
		antinode = self.Brother_Node(node)
		TrueReflection = set(self.NodeCluster(antinode))
		ExpectedReflection = set(self.ExpectedClusReflection(node))
		
		if TrueReflection == ExpectedReflection:
			return
		else:
			badNeighbours = TrueReflection - ExpectedReflection
			missingNeighbours = ExpectedReflection - TrueReflection
			print(badNeighbours)
			print(missingNeighbours)
			for v in badNeighbours:
				self.remove_edge(self.edge(antinode, v))
			for w in missingNeighbours:
				self.add_edge(antinode, w, add_missing=True)
	
	def ReflectSymVal(self, node):
		""" reflect symmetry value of node """
		sym = self.vp.Symmetry[node]
		if sym > 0:
			if sym > 1:
				sym -= 1
			else:
				sym += 1
		return sym		

	def Reflect(self, v):
		""" make a reflection of a given node, including neighbours and properties """
		# Check for node reflection
		pairid = self.vp.pairing[v]
		pair = gt.find_vertex(self, self.vp.pairing, pairid)
		if len(pair) == 2: 	# If their is a reflection of the given node,
							# fix any differences in topology
			self.FixEdgeSymTopo(v)
			
		else:	# generate new reflection of given node
			v2 = self.add_vertex()
			sym = self.ReflectSymVal(v)
			self.vp.Symmetry[v2] = sym
			self.vp.id[v2] = sym + pairid
			for attr in self.attr:
				attr[v2] = attr[v]
			for neighb in self.ReflectNodes(self.Neighbours(v)):
				self.add_edge(v2, neighb)

	def FixNodeSymmetry(self, nodes):
		""" Fix symmetry of relfection of given nodes"""
		# Reflect assymetries on left side onto right side
		for v in nodes:
			self.Reflect(v)

		# remove assymetries on right side
		righties = gt.find_vertex(self, self.vp.Symmetry, 2)
		for v in righties: 
			if len(gt.find_vertex(self, self.vp.pairing, self.vp.pairing[v])) != 2:
				self.remove_vertex(v)

		return

	def FixEdgeSymmetry(self, nodes):
		""" fix symmetry value and topology for 
		all nodes"""
		for v in nodes:
			self.FixEdgeSymTopo(v)
	
	def FixMorphSymmetry(self):
		""" fix symmetry value and topology for 
		all nodes, choosing left nodes as truth"""
		lefties = gt.find_vertex(self, self.vp.Symmetry, 1)
		self.FixNodeSymmetry(lefties)
		#self.FixEdgeSymmetry(lefties)

	def add_node_pair(self):
		""" add pair of nodes """
		v = self.add_vertex()
		u = self.add_vertex()
		pairid = int(u)*100
		self.vp.Symmetry[v] = 1
		self.vp.Symmetry[u] = 2
		self.vp.pairing[v] = pairid
		self.vp.pairing[u] = pairid
		self.vp.id[v] = pairid + 1
		self.vp.id[u] = pairid + 2
		self.vp.brother[u] = v
		self.vp.brother[v] = u
		return v, u

	def add_node(self):
		""" add middle node """
		v = self.add_vertex()
		pairid = int(v)*100
		self.vp.Symmetry[v] = 0
		self.vp.pairing[v] = pairid
		self.vp.id[v] = pairid + 0
		self.vp.brother[v] = v

		return v

	def add_edge_pair(self, u, v):
		av = self.Brother_Node(v)
		au = self.Brother_Node(u)
		self.add_edge_list([(av, au), (u,v)])

	def remove_node_pair(self, v):
		""" symmetrically remove pair of nodes """
		vs = gt.find_vertex(self, self.vp.pairing, self.vp.pairing[v])
		self.remove_vertex(vs)

	def remove_nodes_from(self, nodes):
		""" remove listed nodes from g """
		for v in nodes:
			self.remove_node_pair(v)

	def remove_edge_pair(self, u, v):
		""" symmetrically remove edge pair """
		u2 = self.Brother_Node(u)
		v2 = self.Brother_Node(v)
		self.remove_edge(self.edge(u,v))
		try:
			self.remove_edge(self.edge(u2,v2))
		except:
			pass

	def remove_edges_from(self, edges):
		""" remove edges from edgelist """
		for e in edges:
			self.remove_edge_pair(e[0], e[1])
		return

	def add_edges_from(self, edges):
		""" add edges from edgelist """
		for e in edges:
			self.add_edge_pair(e[0], e[1])

	def choose_ran_node(self, exclude = []):
		""" choose random node that is not of a list of excluded nodes """
		v = exclude[0]

		while v in exclude:
			v = np.random.choice(self.get_vertices())

		return v

	def add_n_edges_to_node(self, v, n=1):
		""" adds n random edges to a given node """
		degree = self.degree(v)
		target = degree + n

		if target > self.num_vertices()-1:
			return 0
		
		# Reflect all actions
		a1 = set(self.Neighbours(v))
		
		
		while self.degree(v) < target:
			# identify second degree neighbours to connect with
			u = self.choose_ran_node(exclude = [v] + list(a1))
			# add edge between target node and 2nd degree neighbour, and reflect
			self.add_edge_pair(v, u)
		
		return 0

	def node_merge(self, u, v):
		""" Merge 2 node pairs """
		neighbs = self.Neighbours(v)
		us = [(u, n) for n in neighbs if n != u]
		self.add_edges_from(us)
		self.remove_node_pair(v)

	def New_Node_Name(self, reflect=False, middle=False):
		""" Create new random node (50% change middle or paired) """
		indices = []
		if reflect:
			indices += list(self.add_node_pair())
		elif middle:
			indices.append(self.add_node())
		elif np.random.randint(0,2):
			indices += list(self.add_node_pair())
		else:
			indices.append(self.add_node())

		return indices

class Placoderm(biMorph):
	""" Adds mutator functions relevant to placoderms """
	def __init__(self, directed=False, prune=False, vorder=None):
		biMorph.__init__(self, directed=False, prune=False, vorder=None)

	def mutate(self, Node, move):
		""" Choose and perform move """
		Move_dict = {
			0:self.char_grows,
			1:self.char_shrinks,
			2:self.char_moves,
			3:self.char_gain,
			4:self.char_loss,
			5:self.char_merge,
			6:self.char_split,
			7:self.char_expansion,
			8:self.char_squeein,
			9:self.char_squeeout
		}
		Move_dict[move](Node)
		return 0

	def mutator(self, move = None, Node = None):
		""" Make random move """
		movement = None
		try_again = True

		while try_again == True:
			#Copy of graph:
			G = self.copy()
			G.__class__ = Placoderm
			
			#Variables and Constraints
			G_Size = G.num_vertices()

			if Node == None:
				#neighbs = nx.single_source_shortest_path_length(G, 'body')
				#neighbs = [k for k in neighbs if  k[0] != "L"]
				Node = np.random.choice(G.get_vertices())

			if move == None:
				moves = list(range(10))
				if G_Size <= 2:              # G can only gain or split
					moves = [0, 1, 2,		3, 6, 7]
				elif G.degree(Node) == 0:      # Node must grow or be lost
					print("boop")
					moves = [0, 1, 2, 		4, 5, 8]
				move = np.random.choice(moves)
				
			G.mutate(Node, move)
			#neighbs = nx.single_source_shortest_path_length(G, np.random.choice(list(G.nodes())))
			#neighbs = [k for k in neighbs]
			
			# Ensure still attached to body
			if 1==1:
				break
			else:
				Node = None
				move = None

		# remove selfloops and parallels
		gt.remove_self_loops(G)  
		gt.remove_parallel_edges(G)

		# remove solitary nodes
		remove = []
		degrees = G.degrees()
		for i in degrees:
			if degrees[i] == 0:
				remove.append(i)
		G.remove_nodes_from(remove)

		return G, move
	
	def char_grows(self, node = None):
		""" Plate grows, so node gains edges """
		if node == None:
			node = np.random.choice(self.get_vertices())

		self.add_n_edges_to_node(v=node)
	
		return 0
	
	def char_shrinks(self, node = None):
		""" Plate shrinks, so loses edges. If leaf node remove """
		if node == None:
			node = np.random.choice(self.get_vertices())
			
		v = np.random.choice(self.Neighbours(node))
		
		#remove collected edges
		self.remove_edge_pair(node, v)
	
	def char_moves(self, node = None):
		""" Plate moves so edges are replaced """
		if node == None:
			node = np.random.choice(self.get_vertices())
		
		self.char_grows(node)
		self.char_shrinks(node)
	
		return 0
	
	def char_gain(self, Node = None):
		""" New plate emerges, new node with mean number edges """
		node1 = self.New_Node_Name()
		self.add_edge_pair(node1[0], Node)
		
		if self.num_vertices() > 2:
			n = self.num_vertices()
			while n >= self.num_vertices():
				# Add typical number of edges
				n = self.ran_degree() 
			n -= 1
			
			if self.vp.Symmetry == 0:
				n /= 2
				n = int(n)
			# reflected nodes
			self.add_n_edges_to_node(node1[0], n=n)
		
		# Correct symmetry label
		self.FixNodeSymmetryValue(node1[0])
		try:
			self.FixNodeSymmetryValue(node1[1])
		except:
			pass
		return 0
	
	def char_loss(self, Node = None):
		""" plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
		
		if not Node:
			Node = np.random.choice(self.get_vertices())
		degree = self.degree(Node)
	
		if degree > 1:
			
			neighbs = np.random.choice(list(self.Neighbours(Node)), \
				np.random.randint(0, self.num_vertices()))

			self.add_edges_from(
				it.product(
					neighbs,
					neighbs
				)
			)

		self.remove_node_pair(Node)
	
		return 0
	
	def char_merge(self, u = None):
		""" 2 plates become one, 2 adjacent nodes become the same node, union of edges """
		
		if u == None:
			u = np.random.choice(self.get_vertices())
		
		neighbs = list(self.Neighbours(u))
		v = str(np.random.choice(neighbs))

		self.node_merge(u,v)
				
		return 0
	
	def char_split(self, u = None):
		""" 0ne plate becomes 2, half the instances of a node 
		are replaced with new node that will be adjacent to old node """

		print(self.vp.id[u])
		uneighbours = self.Neighbours(u)
		ud = len(uneighbours)
		
		usym = self.vp.Symmetry[u]
		if usym:
			v = self.New_Node_Name(reflect=True)[0]
		else:
			v = self.New_Node_Name()[0]
		vsym = self.vp.Symmetry[v]
	
		if ud > 1: #Real split
			if not usym and vsym: # 
				neighbs = self.Neighbours(u)
				self.remove_node_pair(u)
				vpair = self.Brother_Node(v)
				self.add_edge_pair(v, vpair)
				v1s = [(v, n) for n in neighbs \
					if self.vp.Symmetry[n] == vsym or self.vp.Symmetry[n] ==0]
				self.add_edges_from(v1s)

			else: # mid to mods or refs to refs
				# Collect 50% of adjacencies
				neighbs = np.random.choice( uneighbours, size=np.random.randint(0, ud), replace=False )
				
				if not vsym:	# Make sure exactly one middle neighbour
								# if middle node being split into middle 
								# nodes
					midss = set(gt.find_vertex(self, self.vp.Symmetry, 0))
					mids = midss & set(neighbs)
					if len(mids) > 1:
						mid = np.random.choice(list(mids))
						mids -= {mid}
						neighbs = set(neighbs) - mids
					elif len(mids) < 1 and len(neighbs)>0:
						pos_mids = midss & set(uneighbours)
						neighbs[0] = np.random.choice(list(pos_mids))

				for a in neighbs:
					try:
					# Remove adjacencies from old node0
						self.remove_edge_pair(u, a)
					except:
						pass
					# Give adjacencies to new node
					self.add_edge_pair(v, a)

				#Add edge between split nodes
				self.add_edge_pair(u, v)
		
		else:	#Pseudo split, equivalent to char_gain
			self.add_edge_pair(u,v)
		
		return 0
	
	def char_expansion(self, v = None):
		""" Where 2 characters meet, add new character adjacent to them
		and their coneighbours """
		# choose neighbour:
		u = np.random.choice(self.Neighbours(v))
		self.remove_edge_pair(u,v)
		if self.vp.Symmetry[u] + self.vp.Symmetry[v] != 0:
			w = self.New_Node_Name(reflect=True, middle=True)[0]
		else:
			w = self.New_Node_Name(middle=True, reflect=False)[0]

		cneighbs = self.coneighbours([u, v])
		self.add_edge_pair(w, u)
		self.add_edge_pair(w, v)
		pop = len(cneighbs)

		if pop > 0 and np.random.randint(0,3): 
			for i in cneighbs:
				self.add_edge_pair(w, i)

		return 0
	
	def char_squeein(self, node=None):
		""" similar to expansion but moves an existing node in, instead of creating a new one """
		sym = self.vp.Symmetry[node]
		if sym > 0: #L
			S = set(gt.find_vertex(self, self.vp.Symmetry, sym))
			neighbs = set(self.Neighbours(node))
			neighbs = S & neighbs
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(self.Neighbours(neighb1))
			neighbs2 = S & neighbs2
			neighbs2 -= {node}
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))

		else: #M
			neighbs = set(self.Neighbours(node))
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(self.Neighbours(neighb1)) & neighbs
			neighbs2 -= {node}
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))
		
		co_neighbs = set(self.Neighbours(neighb1)) & set(self.Neighbours(neighb2))
		self.add_edge_pair(node, neighb1)
		self.add_edge_pair(node, neighb2)
		self.remove_edge_pair(neighb1, neighb2)


		if co_neighbs:
			co_neighb = np.random.choice(list(co_neighbs))
			self.add_edge_pair(node, co_neighb)

		return 0
	
	def char_squeeout(self, node=None): 
		""" reverse of squeein """

		neighbs = self.Neighbours(node)
		np.random.shuffle(neighbs)
		for neighb1 in neighbs:
			cneighbs = list(self.coneighbours([neighb1, node]))
			if len(cneighbs) > 1:
				self.remove_edge_pair(node, neighb1)
				self.add_edge_pair(cneighbs[0], cneighbs[1])
				return 0
		
		return 0

## Functions ##
def Profiler(G):
	for i in range(1597):
		neighbs = G.mutator(move=np.random.choice(range(9)))

## Main ##
def main(argv):
	""" Test mutator """	
	G1 = Placoderm.From_Dir(argv[1])
	G1,a = G1.mutator(move=np.random.choice(range(9)))

	#cProfile.runctx("Profiler(G1)", {'G1':G1, 'Profiler':Profiler}, {})
	G1.draw()
	G1.list_properties()
	print(G1.vertex_properties['Symmetry'])
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)