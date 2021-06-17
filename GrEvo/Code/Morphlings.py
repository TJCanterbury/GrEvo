#!/usr/bin/env python3

""" The Morphling and daughter classes for evolving from one species to another """

__appname__ = 'Morphlings.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import itertools as it
import numpy as np
import networkx as nx
import gressure as gr
import matplotlib.pyplot as plt
import csv
import re

## Classes ##
class Morphling(nx.Graph):
	""" Adds useful functions to the nx.Graph class """
	def __init__(self, data=None, **attr):
		nx.Graph.__init__(self, data, **attr)
	# General graph functions
	def num_edges(graph):
		""" Returns from a normal distribution an integer for the
		number of edges to add """
		my_degrees = graph.degree()
		degree_values = [v for k, v in my_degrees]
		sum_G = sum(degree_values)

		mean = (sum_G / graph.number_of_nodes())
		randomInts = np.random.normal(loc=mean, size=1).astype(int)
		
		while randomInts < 0:
			randomInts = np.random.normal(loc=mean, size=1).astype(int)
		
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

		if one == True:
			node = gr.Reflect([node])
			node = node[0]

		else:
			node = gr.Reflect(list(node))


		return node

	def add_n_edges_to_node(graph, n, node):
		""" adds n random edges to a given node """
		
		target = graph.degree(node) + n
		c1 = list()
		c2 = list()

		if target > graph.number_of_nodes()-1:
			return graph, c1, c2
		node2 = graph.reflect_n(node)
		
		# Reflect all actions
		a1 = set(graph.neighbors(node))
		a2 = set(graph.neighbors(node2))
		
		while graph.degree(node) < target:
			# identify second degree neighbours to connect with
			neighbs = nx.single_source_shortest_path_length(graph, node, cutoff=2)
			second_neighbs = [k for k in neighbs if neighbs[k] == 2]

			if not second_neighbs:
				v = graph.choose_ran_node(exclude = [node, graph.reflect_n(node)])
			else:
				# pick one
				v = np.random.choice(second_neighbs)

			# add edge between target node and 2nd degree neighbour, and reflect
			graph.add_edge(node, str(v))
			graph.add_edge(node2, graph.reflect_n(str(v)))
		
		b1 = set(graph.neighbors(node))
		b2 = set(graph.neighbors(node2))
		
		c1 = list(b1 - a1)
		c2 = list(b2 - a2)
		
		return graph, c1, c2

	def add_ran_node(graph, reflect=False, middle=False):
		""" Add node """
		node = str(graph.number_of_nodes() + 1)
		while graph.has_node(node):
			node = str(int(node) + 1)

		if (np.random.randint(0,2) or reflect) and not middle:
			new_node = "L" + node
			while graph.has_node(new_node):
				new_node = "L" + str(int(new_node[1:]) + 1)
			graph.add_node(new_node)
			node = new_node

		else:
			graph.add_node(node)

		return graph, node

	def choose_ran_node(graph, exclude = ["body"]):
		""" choose random node that is not of a list of excluded nodes """
		v = exclude[0]

		while v in exclude:
			v = np.random.choice(graph.nodes())

		return v

	def LCR(G1, G2, stuck=False, node="body"):
		""" Find the largest common radius of nodes about the body, 
		to exclude these from morphing """
		
		# Initiate variables
		radius = 0
		neighbs_1 = nx.single_source_shortest_path_length(G1, "body")
		max_radius_1 = max(neighbs_1.values())
		
		if nx.is_isomorphic(G1, G2): # If solution found return score of 0
			return 0, radius

		while radius <= max_radius_1: # Else, start by matching ego_graph about the body
			A = nx.ego_graph(G1, node, radius=radius, undirected= True)
			B = nx.ego_graph(G2, node, radius=radius, undirected= True)
			
			if nx.is_isomorphic(A, B):
				save_rad = radius
				radius += 1
			else:
				break

		return save_rad

	@classmethod
	def from_edgelist(cls, filex):
		""" derives Morphling from edgelist format """
		with open(filex) as infile:
			csv_reader = csv.reader(infile, delimiter=' ')
			return cls(csv_reader)

	def __str__(self, file="display.png"):
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
	
	def perturber(graph, move = None, Node = None, stuck = False):
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
				neighbs = [k for k in neighbs if neighbs[k] >= 1 and k[0] != "L" ]
				Node = np.random.choice(neighbs)

			if move == None:
				moves = list(range(0, 10))
				if G_Size <= 2:              # G can only gain or split
					moves = [3, 6, 7]
				elif G.degree(Node) == 0:      # Node must grow or be lost
					moves = [4, 5, 8]
				if True:
					moves.extend( [0, 1, 2] )
				move = np.random.choice(moves)

			if move == 0: # grows
				G, movement = G.char_grows(Node)
			if move == 1: # shrinks
				G, movement = G.char_shrinks(Node)
			if move == 2: # moves
				G, movement = G.char_moves(Node)
			if move == 3: # gain
				G, movement = G.char_gain(Node)
			if move == 4: # loss
				G, movement = G.char_loss(Node)
			if move == 5: # merge
				G, movement = G.char_merge(Node)
			if move == 6: # split
				G, movement = G.char_split(Node)
			if move == 7: # Graph expansion (replace edge with node)
				G, movement = G.char_expansion(Node)
			if move == 8: # squeeze node between 2 others
				G, movement = G.char_squeein(Node)
			if move == 9: # move node out from between 2 nodes
				G, movement = G.char_squeeout(Node)

			neighbs = nx.single_source_shortest_path_length(G, 'body')
			neighbs = [k for k in neighbs]
			# Ensure still attached to body
			if len(neighbs) > 1 and movement:
				try_again = False
			else:
				G = graph.copy()
				Node = None
				move = None

		# remove selfloops
		G.remove_edges_from(nx.selfloop_edges(G))   

		# remove solitary nodes
		remove = [node for node,degree in dict(G.degree()).items() if degree == 0]
		G.remove_nodes_from(remove)
		
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
		""" plate grows, so node gains edges """
		if node == None:
			node = np.random.choice(graph.nodes())

		node2 = graph.reflect_n(node)
	
		graph, b1, b2 = graph.add_n_edges_to_node(1, node)
	
		edges =  [node + " - " + ', '.join(b1), node2 + " - " +  ', '.join(b2)]
		movement =  "char_grows- Added edges: " + edges[0] + "; " + edges[1]
	
		return graph, movement
	
	def char_shrinks(graph, node = None):
		"""plate shrinks, so loses edges. If leaf node remove """
		if node == None:
			node = np.random.choice(graph.nodes())
		
		# Pick neighbours of node, symmetrically
		nodes = [node, graph.reflect_n(node)]
	
		v = "body"
		while v == 'body':
			v = np.random.choice(list(graph.neighbors(node)))
		vs = [v, graph.reflect_n(v)]
		edges = [(nodes[0], vs[0]), (nodes[1], vs[1])]
		
		#remove collected edges
		graph.remove_edges_from(edges)
		edges = [' - '.join(list(elem)) for elem in edges]
		movement =  "char_shrinks- removed edges: " + edges[0] + "; " + edges[1]
	
		return graph, movement
	
	def char_moves(graph, node = None):
		""" plate moves so edges are replaced """
		if node == None:
			node = str(np.random.choice(graph.nodes()))
		
		graph, m1 = graph.char_grows(node)
		graph, m2 = graph.char_shrinks(node)
	
		movement = "char_moves- node: " + node + \
			". Edge change: " + m2 + " ---> " + m1
	
		return graph, movement
	
	def char_gain(graph, Node = None):
		""" New plate emerges, new node with mean edges of 3 """
		graph, node1 = graph.add_ran_node()
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
			graph, b1, b2 = graph.add_n_edges_to_node(n, node1)
		
		if graph.LR_Ratio(node1):
			new_node1 = graph.LR_Ratio(node1) + graph.add_ran_node(reflect=True)[1][1:]
			new_node2 = graph.reflect_n(new_node1)
			
			mapping = {node1:new_node1, node2:new_node2}
			
			graph = nx.relabel_nodes(graph, mapping)
			node1 = new_node1
			node2 = new_node2

		if node1 != node2:
			movement = "char_gain- Node_1: " + node1 + " new edge(s): " + Node + ", " +  ', '.join(b1) + \
				". Node_2: " + node2 + " new edge(s): " + graph.reflect_n(Node) + ", " + \
				', '.join(b2) + "."
			
		else:
			movement = "char_gain- Node_1: " + node1 + " new edge(s): " + \
				', '.join(b1) + ", " + Node + ", " + graph.reflect_n(Node) + "."
		
		return graph, movement
	
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
		
		movement = "char_loss- node: " + Node + ", " + graph.reflect_n(Node) + "."
	
		return graph, movement
	
	def char_merge(graph, u = None):
		""" 2 plates become one, 2 adjacent nodes become the same node, union of edges """
		
		if u == None:
			u = str(np.random.choice(graph.nodes()))
		
		v = 'body'
		while v == 'body':
			v = str(np.random.choice(list(graph.neighbors(u))))
		
		if u == graph.reflect_n(v):
			node = u[1:]
			graph.add_node(node)
			graph = nx.contracted_nodes(graph, node, u)
			graph = nx.contracted_nodes(graph, node, v)
			
			movement = "char_merge- Eater node: " + node + "." + \
				" eaten node: " + u + ", " + v + "."
	
		elif not graph.is_L_or_R(u) and not graph.is_L_or_R(v): # 2 middle nodes merge
			graph = nx.contracted_nodes(graph, u, v)
			
			movement = "char_merge- Eater node: " + u + "." + " eaten node: " + v + "."
		
		elif not graph.is_L_or_R(u) and graph.is_L_or_R(v): # a middle node merges with 2 side nodes
			graph = nx.contracted_nodes(graph, u, v)
			graph = nx.contracted_nodes(graph, u, graph.reflect_n(v))
			
			movement = "char_merge- Eater node: " + u + "." + " eaten nodes: " + v + ', ' +\
				 graph.reflect_n(v) + "."
	
		elif graph.is_L_or_R(u) and not graph.is_L_or_R(v): # 2 side nodes merge with a middle node
			graph = nx.contracted_nodes(graph, v, u)
			graph = nx.contracted_nodes(graph, v, graph.reflect_n(u))
			
			movement = "char_merge- Eater node: " + v + ". " + \
				"eaten nodes: " + u + ', ' + graph.reflect_n(u) + "."
	
		elif graph.is_L_or_R(u) and graph.is_L_or_R(v): # 2 side nodes merge, 
											# as do their relfections
			graph = nx.contracted_nodes(graph, u, v)
			graph = nx.contracted_nodes(graph, graph.reflect_n(u), graph.reflect_n(v))
		
			movement = \
				"char_merge- Eater nodes: " + u + ", " + graph.reflect_n(u) + \
				". Eaten nodes: " + v + ", " + graph.reflect_n(v) + "."
				
		return graph, movement
	
	def char_split(graph, u1 = None):
		""" 0ne plate becomes 2, half the instances of a node 
		are replaced with new node that will be adjacent to old node """
		if u1 == None:
			u1 = np.random.choice(graph.nodes())
		
		u2 = graph.reflect_n(u1)
		graph, node = graph.add_ran_node()
		v1 = node
		v2 = graph.reflect_n(node)
		ud = graph.degree(u1)
	
		if ud > 1:
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
				
		if graph.is_L_or_R(v1) and not graph.is_L_or_R(u1):
			graph.add_edge(v1, v2)
			graph.remove_node(u1)
		else:
			# Add edge between split nodes
			graph.add_edge(u1, v1)
			graph.add_edge(u2, v2)
		
		movement = "char_split- from node(s): " + u1 + ', ' + u2 + \
			" new node(s): " + v1 + ", " + v2 + "."
	
		return graph, movement
	
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
			graph, v1 = graph.add_ran_node(reflect=True)
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

			movement = "char_expansion- from edge(s): " + u1 + "-" + n1 + ', ' + u2 +  "-" + n2 + \
				" new node(s): " + v1 + ", " + v2 + "."
		
		else:
			# New nodes:
			graph, v1 = graph.add_ran_node(middle=True)
	
			graph.add_edge(u1, v1)
			graph.add_edge(n1, v1)
	
			movement = "char_expansion- from edge(s): " + u1 + "-" + n1 + ', ' + u2 +  "-" + n2 + \
				" new node(s): " + v1 + "."
	
		# Remove old edges
		try:
			graph.remove_edge(u1, n1)
		except:
			pass
		try:
			graph.remove_edge(u2, n2)
		except:
			pass
		
		return graph, movement
	
	def char_squeein(graph, node=None):
		""" similar to expansion but moves an existing node in, instead of creating a new one """

		if graph.is_L_or_R(node) == 2: #R
			neighbs = set(graph.neighbors(node))
			neighbs -= graph.filter_nodes("L")
			if not neighbs:
				return graph, None
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1))
			neighbs2 -= graph.filter_nodes("L")
			neighbs2 -= set(node)
			if not neighbs2:
				return graph, None
			neighb2 = np.random.choice(list(neighbs2))

		if graph.is_L_or_R(node) == 1: #L
			neighbs = set(graph.neighbors(node))
			neighbs -= graph.filter_nodes("R")
			if not neighbs:
				return graph, None
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1))
			neighbs2 -= graph.filter_nodes("R")
			neighbs2 -= set(node)
			if not neighbs2:
				return graph, None
			neighb2 = np.random.choice(list(neighbs2))

		if not graph.is_L_or_R(node): #M
			neighbs = set(graph.neighbors(node))
			if not neighbs:
				return graph, None
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(graph.neighbors(neighb1)) & neighbs
			neighbs2 -= set(node)
			if not neighbs2:
				return graph, None
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
		
		movement = "char_squeein- node: " + node + " moved between nodes: " + neighb1 + "-" + neighb2 + " and node:"  + r_node + " moved between nodes: " + r_neighb1 + "-" + r_neighb2 

		return graph, movement
	
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
			return graph, None
		neighb2 = np.random.choice(list(neighbs))
		r_node = graph.reflect_n(node)
		r_neighb1 = graph.reflect_n(neighb1)
		r_neighb2 = graph.reflect_n(neighb2)

		graph.add_edge(neighb1, neighb2)
		graph.add_edge(r_neighb1, r_neighb2)
		
		movement = "char_squeeout- node: " + node + " squeezed out from between: " + neighb1 + "-" + neighb2 + " and node:"  + r_node + " squeezed out from between: " + r_neighb1 + "-" + r_neighb2 

		return graph, movement

### Business End ###
def main(argv):
	""" Test perturber """
	G1 = Plates.from_edgelist(argv[1])
	morph, move = G1.perturber(Node=argv[2], move=int(argv[3]), radius=1)
	print(move)
	morph.__str__("display_test.png")
	
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)