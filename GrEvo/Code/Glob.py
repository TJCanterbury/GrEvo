#!/usr/bin/env python3

""" Align 2 networks where symmetry has been encoded in node alnmes """

__appalnme__ = 'SymAln.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx
from networkx.algorithms import similarity
from networkx.classes.function import edges, neighbors, nodes
from networkx.convert import from_edgelist
from nxcode import readx
from nxcode import drawx
import numpy as np
import random

## Classes ##
class Morphling(nx.Graph):
	def __init__(graph, data=None, **attr):
		nx.Graph.__init__(graph, data, **attr)
	
	### Basic Morph Functions ### 
	def char_grows(graph, node = None):
		""" plate grows, so node gains edges """
		if node == None:
			node = np.random.choice(graph.nodes())
	
		node2 = reflect_n(node)
	
		graph, b1, b2 = add_n_edges_to_node(graph, 1, node)
	
		edges =  [node + " - " + ', '.join(b1), node2 + " - " +  ', '.join(b2)]
		movement =  "char_grows- Added edges: " + edges[0] + "; " + edges[1]
	
		return graph, movement
	
	def char_shrinks(graph, node = None):
		"""plate shrinks, so loses edges. If leaf node remove """
		if node == None:
			node = np.random.choice(graph.nodes())
		
		# Pick neighbours of node, symmetrically
		nodes = [node, reflect_n(node)]
	
		v = "body"
		while v == 'body':
			v = np.random.choice(list(graph.neighbors(node)))
		vs = [v, reflect_n(v)]
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
		graph, node1 = add_ran_node(graph)
		node2 = reflect_n(node1)
		graph.add_node(node1)
		graph.add_node(node2)
		graph.add_edge(node1, Node)
		graph.add_edge(node2, reflect_n(Node))
		
		if graph.number_of_nodes() > 2:
			n = graph.number_of_nodes()
			while n >= graph.number_of_nodes():
				# Add typical number of edges
				n = num_edges(graph)
			
			# reflected nodes
			graph, b1, b2 = add_n_edges_to_node(graph, n, node1)
			
			Node += ", " + ', '.join(b1) + ". Node_2: " + node2 +\
				 " new edge(s): " + reflect_n(Node) + ", " + ', '.join(b2)
		
		movement = "char_gain- Node_1: " + node1 + " new edge(s): " + Node + "."
		
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
			ref_neighbs = [reflect_n(node) for node in neighbs]
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
		
		if is_L_or_R(Node):
			graph.remove_node(reflect_n(Node))
		graph.remove_node(Node)
		
		movement = "char_loss- node: " + Node + ", " + reflect_n(Node) + "."
	
		return graph, movement
	
	def char_merge(graph, u = None):
		""" 2 plates become one, 2 adjacent nodes become the same node, union of edges """
		
		if u == None:
			u = str(np.random.choice(graph.nodes()))
		
		v = 'body'
		while v == 'body':
			v = str(np.random.choice(list(graph.neighbors(u))))
		
		if u == reflect_n(v):
			node = u[1:]
			graph.add_node(node)
			graph = nx.contracted_nodes(graph, node, u)
			graph = nx.contracted_nodes(graph, node, v)
			
			movement = "char_merge- Eater node: " + node + "." + \
				" eaten node: " + u + ", " + v + "."
	
		elif not is_L_or_R(u) and not is_L_or_R(v): # 2 middle nodes merge
			graph = nx.contracted_nodes(graph, u, v)
			
			movement = "char_merge- Eater node: " + u + "." + " eaten node: " + v + "."
		
		elif not is_L_or_R(u) and is_L_or_R(v): # a middle node merges with 2 side nodes
			graph = nx.contracted_nodes(graph, u, v)
			graph = nx.contracted_nodes(graph, u, reflect_n(v))
			
			movement = "char_merge- Eater node: " + u + "." + " eaten nodes: " + v + ', ' +\
				 reflect_n(v) + "."
	
		elif is_L_or_R(u) and not is_L_or_R(v): # 2 side nodes merge with a middle node
			graph = nx.contracted_nodes(graph, v, u)
			graph = nx.contracted_nodes(graph, v, reflect_n(u))
			
			movement = "char_merge- Eater node: " + v + ". " + \
				"eaten nodes: " + u + ', ' + reflect_n(u) + "."
	
		elif is_L_or_R(u) and is_L_or_R(v): # 2 side nodes merge, 
											# as do their relfections
			graph = nx.contracted_nodes(graph, u, v)
			graph = nx.contracted_nodes(graph, reflect_n(u), reflect_n(v))
		
			movement = \
				"char_merge- Eater nodes: " + u + ", " + reflect_n(u) + \
				". Eaten nodes: " + v + ", " + reflect_n(v) + "."
				
		return graph, movement
	
	def char_split(graph, u1 = None):
		""" 0ne plate becomes 2, half the instances of a node 
		are replaced with new node that will be adjacent to old node """
		if u1 == None:
			u1 = np.random.choice(graph.nodes())
		
		u2 = reflect_n(u1)
		graph, node = add_ran_node(graph)
		v1 = node
		v2 = reflect_n(node)
		ud = int(graph.degree(str(u1)))
	
		if ud > 1:
			# Collect 50% of adjacencies
			neighbs = set(np.random.choice(list(graph.neighbors(u1)), size=np.random.randint(1, ud)))
			neighbs = list(neighbs)
			for a in neighbs:
				a2 = reflect_n(a)
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
				
		if is_L_or_R(v1) and not is_L_or_R(u1):
			graph.add_edge(v1, v2)
			graph.remove_node(u1)
		else:
			# Add edge between split nodes
			graph.add_edge(u1, v1)
			graph.add_edge(u2, v2)
		
		movement = "char_split- from node(s): " + u1 + ', ' + u2 + " new node(s): " + v1 + ", " + v2 + "."
	
		return graph, movement
	
	def char_expansion(graph, u1 = None):
		""" Add node to edge """
		u2 = reflect_n(u1)
		# choose neighbour:
		neighbs1 = list(graph.neighbors(u1))
		n1 = np.random.choice(list(neighbs1), size=1)[0]
		neighbs2 = list(graph.neighbors(n1))

		common_neighbs = set(neighbs1).intersection( set(neighbs2) )
		if len(common_neighbs) > 0:
			reflected_common_neighbs = reflect_n(common_neighbs, one=False)

		n2 = reflect_n(n1)
		
		if is_L_or_R(n1) or is_L_or_R(u1):
			# New nodes:
			graph, v1 = add_ran_node(graph, reflect=True)
			v2 = reflect_n(v1)
	
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
			graph, v1 = add_ran_node(graph, middle=True)
	
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
	
	def char_squeein(graph, node=None): #Incomplete
		return None
	
	def char_squeeout(graph, node=None): #Incomplete
		return None
	
	def seg_duplicate(graph, node=None): #Incomplete

		return 0

	def __str__(graph, file="display.png"):
		drawx(graph, file)
		return 0

### Align and Measure Similarity ###
def reflect_n(node):
	""" returns the reflected node label """

	if node[0] == 'L':
		# If on left side add copy on right side
		node = 'R' + node[1:]
	
	elif node[0] == 'R':
		# If on right side add copy on left side
		node = 'L' + node[1:] 
	
	return node

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

def add_ran_node(graph, reflect=False, middle=False):
	""" Add node """
	node = str(graph.number_of_nodes() + 1)
	while graph.has_node(node):
		node = str(int(node) + 1)
		
	if (np.random.randint(0,2) or reflect) and not middle:
		node = "L" + node
		graph.add_node(node)
		
		return graph, node

	else:
		graph.add_node(node)

		return graph, node

### Align ###
def S3(G1, G2, aln):
	edges_1 = set(G1.edges())
	edges_2 = set(G2.edges())
	c = 0
	d = 0
	for e in edges_1:
		if G2.has_edge(aln[e[0]], aln[e[1]]):
			c += 1
	for f in edges_2:
		if G1.has_edge(f[0], f[1]):
			d += 1

	score = c / ( len(edges_1) - c + len(edges_2))

	return score

def Node_sim(G1, G2, node1, node2, aln):
	D_1 = G1.degree(node1)
	D_2 = G2.degree(node2)
	N_1 = G1.neighbors(node1)
	N_2 = G2.neighbors(node2)
	Sym_1 = is_L_or_R(node1)
	Sym_2 = is_L_or_R(node2)
	if Sym_1 == Sym_2:
		score = 0.5
	else:
		score = 0
	if not D_1 or not D_2:
		score = 0.5
		if Sym_1 == Sym_2:
			score += 0.25
		return score
	Aligned_neighbs = [ neighb for neighb in N_1 if neighb in aln and aln[neighb] in N_2 ]
	score += len( Aligned_neighbs ) / D_1

	return score

def extend_aln(G1, G2, aln, n1, n2):
	""" attempt to add nodes to outer perimeter of subgraphs that 
	maintain isomorphism between these subgraphs and so building 
	the largest common subgraph """

	# for each of n1, align the node with the node in n2 with the highest similairty
	for node in n1:
		# build dictionary of similairites between node and all available 
		# nodes in the opposing graph
		if aln[node]:
			break
		Similarities = { node2: Node_sim(G1, G2, node, node2, aln) for node2 in n2 }
		if Similarities:
			# Assign alignment of first best node based on similarity scores
			aln[node] = max(Similarities, key=Similarities.get)
			aln[reflect_n(node)] = reflect_n(aln[node])

			# remove aligned nodes from subset that are unaligned
			n1.remove(node)
			try: # if node has a reflection it will be removed too
				n2.remove(aln[node])
			except:
				pass
			try: # if node has a reflection it will be removed too
				n1.remove(reflect_n(node))
			except:
				pass
			try: # if node has a reflection it will be removed too
				n2.remove(reflect_n(aln[node]))
			except:
				pass
	return aln

def ran_align(G1, G2, source=[["body", "body"]]):

	# start alignment with the known node correspondance(s) from source
	aln = dict.fromkeys(G1.nodes())
	for i in source:
		aln[i[0]] = i[1]
	
	# make independent alignment from each node in source
	n1 = nx.single_source_shortest_path_length(G1, source[0][0])
	n2 = nx.single_source_shortest_path_length(G2, source[0][1])

	max_radius = min([max(n1.values()), max(n2.values())])
	radius = 1

	while radius <= max_radius:
		A_neighbs = [k for k in n1 if n1[k] == radius] 
		B_neighbs = [k for k in n2 if n2[k] == radius] 
		random.shuffle(A_neighbs)
		random.shuffle(B_neighbs)
		aln = extend_aln(G1, G2, aln, A_neighbs, B_neighbs)
		radius +=1
	
	n1 = [node for node in G1.nodes() if not aln[node]]
	n2 = [node for node in G2.nodes() if node not in aln.values()]
	aln = extend_aln(G1, G2, aln, n1, n2)

	
	# choose alignment with best symmetric substructure score
	return S3(G1, G2, aln), aln

def Align(G1, G2, repeats=10, source=[["body", "body"]]):
	
	# Make the graphs the same size so that a global alignment can be made
	if G1.number_of_nodes() > G2.number_of_nodes():
		G1, G2 = G2, G1
	while G1.number_of_nodes() < G2.number_of_nodes():
		G1, node = add_ran_node(G1, middle=True)
	
	best_score = 0
	best_aln = None
	
	for i in range(repeats):
		score, aln = ran_align(G1.copy(), G2.copy(), source)
		if score > best_score:
			best_score = score
			best_aln = aln
	
	return best_score, best_aln

### Main ###
def main(argv):
	G1 = readx(argv[1])
	G2 = readx(argv[2])
	G1 = Morphling(G1)
	G2 = Morphling(G2)



	score, aln = Align(G1, G2, repeats=500)
	print(aln)
	print(score)
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)