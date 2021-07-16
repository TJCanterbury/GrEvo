#!/usr/bin/env python3

""" The Morphling and daughter classes for evolving from one species to another """

__appname__ = 'Morphlings.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import csv

## Classes ##
class Morphling(nx.Graph):
	""" Adds useful functions to the nx.Graph class """
	def __init__(self, data=None, **attr):
		nx.Graph.__init__(self, data, **attr)
	# General graph functions
	def num_edges(graph):
		""" Returns from a normal distribution an integer for the
		number of edges to add """
		randomInts = np.random.normal(loc=0, scale=3, size=1).astype(int)
		max_degree = sorted(graph.degree, key=lambda x: x[1], reverse=True)[0][1]
		while randomInts < 0 or randomInts > max_degree:
			randomInts = np.random.normal(loc=0, size=1).astype(int)
		
		return int(randomInts)

	@staticmethod
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

	@staticmethod
	def reflect_n(node, one = True):
		""" returns the reflected node label """

		if one:
			if node[0] == "L":
				return "R" + node[1:]
			elif node[0] == "R":
				return "L" + node[1:]
			else:
				return node

		else:
			nodes = list(node)
			r_nodes = []
			for v in nodes:
				if v[0] == "L":
					r_nodes.append("R" + v[1:])
				elif v[0] == "R":
					r_nodes.append("L" + v[1:])
				else:
					r_nodes.append(v)

			return r_nodes

	def copy_attrs(graph, source, sink):
		""" Give new nodes the parent nodes' attributes """
		attrs = {sink: graph.nodes[source]}
		nx.set_node_attributes(graph, attrs)

	def node_merge(graph, u, v):
		""" Merge 2 nodes """
		neighbs = list(graph.neighbors(v))
		us = [(u, n) for n in neighbs if n != u]
		graph.add_edges_from(us)
		graph.copy_attrs(v, u)
		graph.remove_node(v)

	def nodes_merge(graph, u, lv):
		""" Merge multiple nodes with a given node"""
		neighbs = []
		for i in lv:
			neighbs += list(graph.neighbors(i))
			graph.copy_attrs(i, u)
			graph.remove_node(i)
		us = [(u, n) for n in neighbs if n not in lv]
		graph.add_edges_from(us)

	def add_n_edges_to_node(graph, n, node):
		""" adds n random edges to a given node """
		
		target = graph.degree(node) + n

		if target > graph.number_of_nodes()-1:
			return 0
		node2 = graph.reflect_n(node)
		
		# Reflect all actions
		a1 = set(graph.neighbors(node))
		
		
		while graph.degree(node) < target:
			# identify second degree neighbours to connect with
			v = graph.choose_ran_node(exclude = [node] + list(a1))
			# add edge between target node and 2nd degree neighbour, and reflect
			graph.add_edge(node, str(v))
			graph.add_edge(node2, graph.reflect_n(str(v)))
		
		return 0

	def New_Node_Name(graph, reflect=False, middle=False):
		""" Create new node name """
		node = str(graph.number_of_nodes() + 1)
		while graph.has_node(node):
			node = str(int(node) + 1)

		if (np.random.randint(0,2) or reflect) and not middle:
			new_node = "L" + node
			while graph.has_node(new_node):
				new_node = "L" + str(int(new_node[1:]) + 1)
			node = new_node

		return node

	def choose_ran_node(graph, exclude = ["body"]):
		""" choose random node that is not of a list of excluded nodes """
		v = exclude[0]

		while v in exclude:
			v = np.random.choice(graph.nodes())

		return v

	def homologs(graph):
		""" return number of unique homologies """
		ua = nx.get_node_attributes(graph, "homolog")
		return set(ua.values())

	@classmethod
	def from_edgelist(cls, filex):
		""" derives Morphling from edgelist format """
		with open(filex) as infile:
			csv_reader = csv.reader(infile, delimiter=' ')
			return cls(csv_reader)
	
	def attr_from_csv(graph, file):
		""" Set attributes from csv file, each row being 'node key value' """
		attrs = pd.read_csv(file, sep=" ", index_col=0)
		attr = attrs.to_dict(orient='index')
		nx.set_node_attributes(graph, attr)
		return 0

	def __str__(self, file="display.jpeg"):
		""" Saves Morphling graph to a given image file """
		d = dict(self.degree)
		plt.figure(figsize=(7,7))
		nx.draw_networkx(self, with_labels = True, node_color="#d9ffb3", \
		edgecolors='black', node_size=[v * 200 for v in d.values()])
		plt.axis('off')
		plt.savefig(file)
		plt.close()
		return 0

class Plates(Morphling):
	""" Morphling moves for a network of tectonic plate like characters, 
	with physical contact between characters giving an edge """
	def __init__(self, data=None, **attr):
		Morphling.__init__(self, data, **attr)

	def mutate(G, Node, move):
		""" Choose and perform move """
		Move_dict = {
			0:G.char_grows,
			1:G.char_shrinks,
			2:G.char_moves,
			3:G.char_gain,
			4:G.char_loss,
			5:G.char_merge,
			6:G.char_split,
			7:G.char_expansion,
			8:G.char_squeein,
			9:G.char_squeeout
		}
		Move_dict[move](Node)
		return 0

	def mutator(graph, move = None, Node = None, stuck = False):
		""" Make random move """
		movement = None
		try_again = True

		while try_again == True:
			#Copy of graph:
			G = graph.copy()

			#Variables and Constraints
			G_Size = G.number_of_nodes()

			if Node == None:
				neighbs = nx.single_source_shortest_path_length(G, 'body')
				neighbs = [k for k in neighbs if  k[0] != "L"]
				Node = np.random.choice(neighbs)
			
			conserve = G.homologs()

			if move == None:
				moves = list(range(0, 10))
				if G_Size <= 2:              # G can only gain or split
					moves = [3, 6, 7]
				elif G.degree(Node) == 0:      # Node must grow or be lost
					moves = [4, 5, 8]
				if True:
					moves.extend( [0, 1, 2] )
				move = np.random.choice(moves)
				
			G.mutate(Node, move)
			
			neighbs = nx.single_source_shortest_path_length(G, np.random.choice(list(G.nodes())))
			neighbs = [k for k in neighbs]
			
			# Ensure still attached to body
			if len(neighbs) == G.number_of_nodes() and conserve == G.homologs() and G.has_node("body"):
				try_again = False
			else:
				Node = None
				move = None

		# remove selfloops
		G.remove_edges_from(nx.selfloop_edges(G))   

		# remove solitary nodes
		remove = [node for node,degree in dict(G.degree()).items() if degree == 0]
		G.remove_nodes_from(remove)
		movement = move
		return G, movement

	### Basic Morph Functions ### 
	def filter_nodes(graph, Initial = None):
		""" returns set of nodes that start with the given letter (L or R) or else, all middle nodes """
		tb_filtered = list(graph.nodes())
		if Initial:
			filtered = [node for node in tb_filtered if node[0] == Initial]
		
		else:
			filtered = [node for node in tb_filtered if node[0] != "R" and node[0] != "L"]

		return set(filtered)

	def LR_Ratio(graph, Node):
		""" find the ratio of left to right neighbours for a given node """
		Neighbours = set(graph.neighbors(Node)) - set(graph.reflect_n(Node))
		L = sum(node[0] == "L" for node in Neighbours)
		R = sum(node[0] == "R" for node in Neighbours)

		if L > R:
			return "L"

		if L < R:
			return "R"
		
		else:
			return None

	def char_grows(graph, node = None):
		""" Plate grows, so node gains edges """
		if node == None:
			node = np.random.choice(graph.nodes())
	
		graph.add_n_edges_to_node(1, node)
	
		return 0
	
	def char_shrinks(graph, node = None):
		""" Plate shrinks, so loses edges. If leaf node remove """
		if node == None:
			node = np.random.choice(graph.nodes())
		
		# Pick neighbours of node, symmetrically
		nodes = [node, graph.reflect_n(node)]
	
		v = np.random.choice(list(graph.neighbors(node)))
		vs = [v, graph.reflect_n(v)]
		edges = [(nodes[0], vs[0]), (nodes[1], vs[1])]
		
		#remove collected edges
		graph.remove_edges_from(edges)
	
		return 0
	
	def char_moves(graph, node = None):
		""" Plate moves so edges are replaced """
		if node == None:
			node = str(np.random.choice(graph.nodes()))
		
		graph.char_grows(node)
		graph.char_shrinks(node)
	
		return 0
	
	def char_gain(graph, Node = None):
		""" New plate emerges, new node with mean number edges """
		node1 = graph.New_Node_Name()
		node2 = graph.reflect_n(node1)

		graph.add_edge(node1, Node)
		graph.add_edge(node2, graph.reflect_n(Node))
		
		if graph.number_of_nodes() > 2:
			n = graph.number_of_nodes()
			while n >= graph.number_of_nodes():
				# Add typical number of edges
				n = graph.num_edges() 
			n -= 1
			
			if not graph.is_L_or_R(node1):
				n /= 2
				n = int(n)
			# reflected nodes
			graph.add_n_edges_to_node(n, node1)
		
		# Correct symmetry label
		if graph.LR_Ratio(node1):
			new_node1 = graph.LR_Ratio(node1) + graph.New_Node_Name(reflect=True)[1:]
			new_node2 = graph.reflect_n(new_node1)
			
			mapping = {node1:new_node1, node2:new_node2}
			
			graph = nx.relabel_nodes(graph, mapping)
		
		return 0
	
	def char_loss(graph, Node = None):
		""" plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
		
		neighbs = nx.single_source_shortest_path_length(graph, 'body')
		neighbs = [k for k in neighbs if neighbs[k] > 1]
		if neighbs:
			Node = np.random.choice(neighbs)
		else:
			Node = np.random.choice(graph.nodes())
		degree = graph.degree(Node)
	
		if degree > 1:
			
			neighbs = np.random.choice(list(graph.neighbors(Node)), \
				np.random.randint(0, graph.number_of_nodes()))
			ref_neighbs = [graph.reflect_n(node) for node in neighbs]
			graph.add_edges_from(
				it.product(
					neighbs,
					neighbs
				)
			)
	
			graph.add_edges_from(
				it.product(
					ref_neighbs,
					ref_neighbs
				)
			)
		
		if graph.is_L_or_R(Node):
			graph.remove_node(graph.reflect_n(Node))
		graph.remove_node(Node)
	
		return 0
	
	def char_merge(graph, u = None):
		""" 2 plates become one, 2 adjacent nodes become the same node, union of edges """
		
		if u == None:
			u = str(np.random.choice(graph.nodes()))
		
		neighbs = list(graph.neighbors(u))
		if graph.is_L_or_R(u):
			neighbs += [graph.reflect_n(u)]
		v = str(np.random.choice(neighbs))
		
		if u == graph.reflect_n(v):
			node = u[1:]
			graph.add_node(node)
			graph.nodes_merge(node, [u, v])
	
		elif not graph.is_L_or_R(u) and not graph.is_L_or_R(v): # 2 middle nodes merge
			graph.node_merge(u, v)
		
		elif not graph.is_L_or_R(u) and graph.is_L_or_R(v): # a middle node merges with 2 side nodes
			graph.node_merge(u, v)
			graph.node_merge(u, graph.reflect_n(v))
	
		elif graph.is_L_or_R(u) and not graph.is_L_or_R(v): # 2 side nodes merge with a middle node
			graph.node_merge(v, u)
			graph.node_merge(v, graph.reflect_n(u))
	
		elif graph.is_L_or_R(u) and graph.is_L_or_R(v): # 2 side nodes merge, 
											# as do their relfections
			graph.node_merge(u, v)
			graph.node_merge(graph.reflect_n(u), graph.reflect_n(v))
				
		return 0
	
	def char_split(graph, u1 = None):
		""" 0ne plate becomes 2, half the instances of a node 
		are replaced with new node that will be adjacent to old node """
		if u1 == None:
			u1 = np.random.choice(graph.nodes())
		
		u2 = graph.reflect_n(u1)
		v1 = graph.New_Node_Name()
		v2 = graph.reflect_n(v1)
		graph.add_node(v1)
		graph.add_node(v2)
		ud = graph.degree(u1)

		# Give new nodes the parent nodes' attributes
		attrs = {v1: graph.nodes[u1], v2: graph.nodes[u2]}
		nx.set_node_attributes(graph, attrs)
	
		if ud > 1:
			if graph.is_L_or_R(v1) and not graph.is_L_or_R(u1):
				neighbs = list(graph.neighbors(u1))
				v1s = [(v1, n) for n in neighbs if n[0] != v2[0]]
				v2s = [(v2, n) for n in neighbs if n[0] != v1[0]]
				graph.add_edges_from(v1s)
				graph.add_edges_from(v2s)
				graph.add_edge(v1, v2)

			
		
			# Collect 50% of adjacencies
			neighbs = np.random.choice( list(graph.neighbors(u1)), size=np.random.randint(0, ud), replace=False )

			for a in neighbs:
				a2 = graph.reflect_n(a)

				# Give adjacencies to new node
				graph.add_edge(v1, a)
				graph.add_edge(v2, a2)
	
				# Remove adjacencies from old node0
				try:
					graph.remove_edge(u1, a)
				except:
					pass
				try:
					graph.remove_edge(u2, a2)
				except:
					pass
				
				#Add edge between split nodes
				graph.add_edge(u1, v1)
				graph.add_edge(u2, v2)
			
			if graph.is_L_or_R(v1) and not graph.is_L_or_R(u1):
				graph.remove_node(u1)
	
		return 0
	
	def char_expansion(graph, u1 = None):
		""" Add node to edge """
		u2 = graph.reflect_n(u1)
		# choose neighbour:
		neighbs1 = list(graph.neighbors(u1))
		n1 = np.random.choice(list(neighbs1), size=1)[0]
		neighbs2 = list(graph.neighbors(n1))

		common_neighbs = set(neighbs1).intersection( set(neighbs2) )
		if len(common_neighbs) > 0:
			reflected_common_neighbs = graph.reflect_n(common_neighbs, one=False)

		n2 = graph.reflect_n(n1)
		
		if graph.is_L_or_R(n1) or graph.is_L_or_R(u1):
			# New nodes:
			v1 = graph.New_Node_Name(reflect=True)
			v2 = graph.reflect_n(v1)
	
			graph.add_edge(u1, v1)
			graph.add_edge(n1, v1)
			#reflect
			graph.add_edge(u2, v2)
			graph.add_edge(n2, v2)

			if len(common_neighbs) > 0:
				for i in common_neighbs:
					graph.add_edge(v1, i)
				for j in reflected_common_neighbs:
					graph.add_edge(v2, j)
		
		else:
			# New nodes:
			v1 = graph.New_Node_Name(middle=True)
	
			graph.add_edge(u1, v1)
			graph.add_edge(n1, v1)
	
		# Remove old edges
		try:
			graph.remove_edge(u1, n1)
		except:
			pass
		try:
			graph.remove_edge(u2, n2)
		except:
			pass
		
		return 0
	
	def char_squeein(graph, node=None):
		""" similar to expansion but moves an existing node in, instead of creating a new one """

		if graph.is_L_or_R(node) == 2: #R
			neighbs = set(graph.neighbors(node))
			neighbs -= graph.filter_nodes("L")
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1))
			neighbs2 -= graph.filter_nodes("L")
			neighbs2 -= set(node)
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))

		if graph.is_L_or_R(node) == 1: #L
			neighbs = set(graph.neighbors(node))
			neighbs -= graph.filter_nodes("R")
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1))
			neighbs2 -= graph.filter_nodes("R")
			neighbs2 -= set(node)
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))

		if not graph.is_L_or_R(node): #M
			neighbs = set(graph.neighbors(node))
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1)) & neighbs
			neighbs2 -= set(node)
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))
		
		co_neighbs = set(graph.neighbors(neighb1)) & set(graph.neighbors(neighb2))
		graph.add_edge(node, neighb1)
		graph.add_edge(node, neighb2)
		graph.remove_edge(neighb1, neighb2)

		#reflect
		r_node = graph.reflect_n(node)
		r_neighb1 = graph.reflect_n(neighb1)
		r_neighb2 = graph.reflect_n(neighb2)
		graph.add_edge(r_node, r_neighb1)
		graph.add_edge(r_node, r_neighb2)

		if co_neighbs:
			co_neighb = np.random.choice(list(co_neighbs))
			graph.add_edge(node, co_neighb)
			graph.add_edge(r_node, graph.reflect_n(co_neighb))
		
		try:
			graph.remove_edge(r_neighb1, r_neighb2)
		except:
			pass

		return 0
	
	def char_squeeout(graph, node=None): 
		""" reverse of squeein """

		neighbs = set(graph.neighbors(node))
		neighb1 = np.random.choice(list(neighbs))
		neighbs -= {neighb1}
		removenodes = []
		for i in neighbs:
			if graph.has_edge(i, neighb1):
				removenodes += i

		neighbs -= set(removenodes)

		if not neighbs:
			return 0
		neighb2 = np.random.choice(list(neighbs))
		r_neighb1 = graph.reflect_n(neighb1)
		r_neighb2 = graph.reflect_n(neighb2)

		graph.add_edge(neighb1, neighb2)
		graph.add_edge(r_neighb1, r_neighb2)
		
		return 0

	def add_symmetry(graph):
		for node in graph.nodes():
			graph.nodes[node]["Symmetry"] = graph.is_L_or_R(node)
		return 0

### Business End ###
def main(argv):
	""" Test mutator """
	G1 = Plates.from_edgelist(argv[1])
	morph, move = G1.mutator( move=int(argv[2]))
	print(move)
	print(morph.nodes(data=True))
	morph.__str__("display_test.jpeg")
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)