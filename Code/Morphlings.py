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
	def ran_degree(self):
		""" Returns from a normal distribution an integer for the
		number of edges to add """
		randomInts = np.random.normal(loc=0, scale=3, size=1).astype(int)
		max_degree = sorted(self.degree, key=lambda x: x[1], reverse=True)[0][1]
		while randomInts < 0 or randomInts > max_degree:
			randomInts = np.random.normal(loc=0, size=1).astype(int)
		
		return int(randomInts)

	def copy_attrs(self, source, sink):
		""" Give new nodes the parent nodes' attributes """
		attrs = {sink: self.nodes[source]}
		nx.set_node_attributes(self, attrs)

	def choose_ran_node(self, exclude = ["body"]):
		""" choose random node that is not of a list of excluded nodes """
		v = exclude[0]

		while v in exclude:
			v = np.random.choice(self.nodes())

		return v

	def homologs(self):
		""" return number of unique homologies """
		ua = nx.get_node_attributes(self, "homolog")
		return set(ua.values())

	def Find_V(self, attr='Sym', attrValue = 0):
		""" Find vertices with given node attr value """
		nodes = [x for x,y in self.nodes(data=True) if y[attr]==attrValue]
		return nodes
	
	def coneighbours(self, nodes):
		""" returns set of common neighbours between 2 nodes """
		cneighbs = set(self._adj[nodes[0]])
		for v in nodes:
			cneighbs = cneighbs & set(self._adj[v])
		return cneighbs

	def draw(self, file="display.jpeg"):
		""" Saves Morphling graph to a given image file """
		d = dict(self.degree)
		plt.figure(figsize=(7,7))
		nx.draw_networkx(self, with_labels = True, node_color="#d9ffb3", \
		edgecolors='black', node_size=[v * 200 for v in d.values()])
		plt.axis('off')
		plt.savefig(file)
		plt.close()
		return 0

class biMorph(Morphling):
	""" Adds useful bilateral symmetry functions to the morphling class """
	def __init__(self, directed=False, prune=False, vorder=None):
		Morphling.__init__(self, directed=False, prune=False, vorder=None)
	
	def Set_Node_Symmetry(self, node):
		""" Assign symmetry VertexProperty values to node """
		sym = int(node) % 10
		self.nodes[node]['Sym'] = sym
		self.nodes[node]['Group'] = int(node) - self.nodes[node]['Sym']
		self.nodes[node]['Bro'] = self.find_brother(node)
	
	def Initialise_Symmetry(self):
		""" Set symmetry for all nodes """
		for v in self.nodes():
			self.Set_Node_Symmetry(v)
	
	def DetectNeighbSymmetry(self, v):
		""" Detect Symmetry based on neighbour symmetry """
		sym = 0
		for w in self._adj[v]:
			if self.nodes[w]['Sym'] == 1:
				sym -= 1
			elif self.nodes[w]['Sym'] == 2:
				sym += 1
		if sym > 0:
			return 2
		if sym < 0: 
			return 1
		return 0
	
	def FixNodeSymmetryValue(self, v):
		""" Set or correct Symmetry of a node """
		self.nodes[v]['Sym'] = self.DetectNeighbSymmetry(v)

	def Brother_V(self, node):
		""" return index of nodes reflection """
		return self.nodes[node]['Bro']

	def find_brother(self, node):
		""" return index of nodes reflection """
		sym = self.nodes[node]['Sym']
		if sym > 0:
			if sym > 1:
				antId = int(node) - 1
			else:
				antId = int(node) + 1
			return antId
		else:
			return node
		
	def ReflectNodes(self, nodes):
		""" return the reflected node names """
		antinodes = []
		for v in nodes:
			antinodes.append(self.Brother_V(v))
		return antinodes

	def add_v_pair(self):
		""" add pair of nodes """
		pairid = np.max(self.nodes())+10
		pairid -= pairid % 10
		v = pairid+1
		u = pairid+2
		self.add_node(v)
		self.add_node(u)
		self.nodes[v]['Sym'] = 1
		self.nodes[u]['Sym'] = 2
		self.nodes[v]['Group'] = pairid
		self.nodes[u]['Group'] = pairid
		self.nodes[u]['Bro'] = v
		self.nodes[v]['Bro'] = u
		return v, u

	def add_v(self):
		""" add middle node """
		v = np.max(self.nodes())+10
		v -= v % 10
		self.add_node(v)
		self.nodes[v]['Sym'] = 0
		self.nodes[v]['Group'] = v
		self.nodes[v]['Bro'] = v
		return v

	def add_e_pair(self, u, v):
		""" add edge symmetrically """
		self.add_edge(u, v)
		av = self.Brother_V(v)
		au = self.Brother_V(u)
		self.add_edge(au, av)

	def remove_v_pair(self, v):
		""" symmetrically remove pair of nodes """
		try:
			self.remove_node(self.Brother_V(v))
		except:
			pass
		try:
			self.remove_node(v)
		except:
			pass

	def remove_v_from(self, nodes):
		""" remove listed nodes from g """
		for v in nodes:
			self.remove_v_pair(v)

	def remove_e_pair(self, u, v):
		""" symmetrically remove edge pair """
		u2 = self.Brother_V(u)
		v2 = self.Brother_V(v)
		self.remove_edge(u,v)
		try:
			self.remove_edge(u2,v2)
		except:
			pass

	def remove_e_from(self, edges):
		""" remove edges from edgelist """
		for e in edges:
			self.remove_e_pair(e[0], e[1])
		return

	def add_e_from(self, edges):
		""" add edges from edgelist """
		for e in edges:
			self.add_e_pair(e[0], e[1])

	def add_n_e_to_v(self, v, n=1):
		""" adds n random edges to a given node """
		degree = self.degree(v)
		target = degree + n

		if target > self.number_of_nodes()-1:
			return 0
		
		# Reflect all actions
		a1 = set(self._adj[v])
		
		
		while self.degree(v) < target:
			# identify second degree neighbours to connect with
			u = self.choose_ran_node(exclude = [v] + list(a1))
			# add edge between target node and 2nd degree neighbour, and reflect
			self.add_e_pair(v, u)
		
		return 0

	def v_merge(self, u, v):
		""" Merge 2 node pairs """
		neighbs = list(self._adj[v])
		us = [(u, n) for n in neighbs if n != u]
		self.add_e_from(us)
		self.remove_v_pair(v)

	def New_v_Name(self, reflect=False, middle=False):
		""" Create new random node (50% change middle or paired) """
		indices = []
		if reflect:
			indices += list(self.add_v_pair())
		elif middle:
			indices.append(self.add_v())
		elif np.random.randint(0,2):
			indices += list(self.add_v_pair())
		else:
			indices.append(self.add_v())

		return indices

class Placoderm(biMorph):
	""" Morphling moves for a network of tectonic plate like characters, 
	with physical contact between characters giving an edge """
	def __init__(self, data=None, **attr):
		biMorph.__init__(self, data, **attr)

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

	def mutator(self, common_homologs, move = None, Node = None, stuck = False):
		""" Make random move """
		movement = None
		try_again = True

		while try_again == True:
			#Copy of graph:
			G = self.copy()

			#Variables and Constraints
			G_Size = G.number_of_nodes()

			if Node == None:
				#neighbs = nx.single_source_shortest_path_length(G, 'body')
				#neighbs = [k for k in neighbs if  k[0] != "L"]
				Node = np.random.choice(list(G.nodes()))

			if move == None:
				moves = list(range(10))
				if G_Size <= 2:              # G can only gain or split
					moves = [0, 1, 2,		3, 6, 7]
				elif G.degree(Node) == 0:      # Node must grow or be lost
					print("boop")
					moves = [0, 1, 2, 		4, 5, 8]
				move = np.random.choice(moves)
				
			G.mutate(Node, move)
			
			# Ensure still attached to body
			if common_homologs == G.homologs():
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

	def char_grows(self, node = None):
		""" Plate grows, so node gains edges """
		self.add_n_e_to_v(node, 1)
	
	def char_shrinks(self, node = None):
		""" Plate shrinks, so loses edges. If leaf node remove """
		v = np.random.choice(list(self._adj[node]))
		self.remove_e_pair(v, node)
	
	def char_moves(self, node = None):
		""" Plate moves so edges are replaced """
		self.char_grows(node)
		self.char_shrinks(node)
	
	def char_gain(self, Node = None):
		""" New plate emerges, new node with mean number edges """
		node1 = self.New_v_Name()[0]
		self.add_e_pair(node1, Node)
		
		if self.number_of_nodes() > 2:
			n = self.number_of_nodes()
			while n >= self.number_of_nodes():
				# Add typical number of edges
				n = self.ran_degree() 
			n -= 1
			
			if self.nodes[node1]['Sym'] == 0:
				n /= 2
				n = int(n)
			# reflected nodes
			self.add_n_e_to_v(node1, n=n)
		
		# Correct symmetry label
		self.FixNodeSymmetryValue(node1)
		try:
			self.FixNodeSymmetryValue(node1)
		except:
			pass
		return 0
	
	def char_loss(self, Node = None):
		""" plate lost, node lost along with edges, if not leaf node the node should be replaced with an edge """
		degree = self.degree(Node)
	
		if degree > 1:
			
			neighbs = np.random.choice(list(self._adj[Node]), \
				np.random.randint(0, self.number_of_nodes()))

			self.add_e_from(
				it.product(
					neighbs,
					neighbs
				)
			)

		self.remove_v_pair(Node)
	
		return 0
	
	def char_merge(self, u = None):
		""" 2 plates become one, 2 adjacent nodes become the same node, union of edges """
		neighbs = list(self._adj[u])
		v = np.random.choice(neighbs)

		self.v_merge(u,v)
	
	def char_split(self, u = None):
		""" 0ne plate becomes 2, half the instances of a node 
		are replaced with new node that will be adjacent to old node """
		uneighbours = self._adj[u]
		ud = len(uneighbours)
		
		usym = self.nodes[u]['Sym']
		if usym:
			v = self.New_v_Name(reflect=True)[0]
		else:
			v = self.New_v_Name()[0]
		vsym = self.nodes[v]['Sym']
	
		if ud > 1: #Real split
			if not usym and vsym: # 
				neighbs = self._adj[u]
				self.remove_v_pair(u)
				vpair = self.Brother_V(v)
				self.add_e_pair(v, vpair)
				v1s = [(v, n) for n in neighbs \
					if self.nodes[n]['Sym'] == vsym or self.nodes[n]['Sym'] ==0]
				self.add_e_from(v1s)

			else: # mid to mods or refs to refs
				# Collect 50% of adjacencies
				neighbs = np.random.choice( list(uneighbours), size=np.random.randint(0, ud), replace=False )
				
				if not vsym:	# Make sure exactly one middle neighbour
								# if middle node being split into middle 
								# nodes
					midss = set(self.Find_V())
					mids = midss & set(neighbs)
					if len(mids) > 1:
						mid = np.random.choice(list(mids))
						mids -= {mid}
						neighbs = set(neighbs) - mids
					elif len(mids) < 1 and len(neighbs)>0:
						pos_mids = midss & set(uneighbours)
						if pos_mids:
							neighbs[0] = np.random.choice(list(pos_mids))

				for a in neighbs:
					try:
					# Remove adjacencies from old node0
						self.remove_e_pair(u, a)
					except:
						pass
					# Give adjacencies to new node
					self.add_e_pair(v, a)

				#Add edge between split nodes
				self.add_e_pair(u, v)
		
		else:	#Pseudo split, equivalent to char_gain
			self.add_e_pair(u,v)
		
		return 0
	
	def char_expansion(self, v = None):
		""" Where 2 characters meet, add new character adjacent to them
		and their coneighbours """
		# choose neighbour:
		u = np.random.choice(list(self._adj[v]))
		self.remove_e_pair(u,v)
		if self.nodes[u]['Sym'] + self.nodes[v]['Sym'] != 0:
			w = self.New_v_Name(reflect=True, middle=True)[0]
		else:
			w = self.New_v_Name(middle=True, reflect=False)[0]

		cneighbs = self.coneighbours([u, v])
		self.add_e_pair(w, u)
		self.add_e_pair(w, v)
		pop = len(cneighbs)

		if pop > 0 and np.random.randint(0,3): 
			for i in cneighbs:
				self.add_e_pair(w, i)

		return 0
	
	def char_squeein(self, node=None):
		""" similar to expansion but moves an existing node in, instead of creating a new one """
		sym = self.nodes[node]['Sym']
		if sym > 0: #L
			S = set(self.Find_V(attrValue=sym))
			neighbs = set(self._adj[node])
			neighbs = S & neighbs
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(self._adj[neighb1])
			neighbs2 = S & neighbs2
			neighbs2 -= {node}
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))

		else: #M
			neighbs = set(self._adj[node])
			if not neighbs:
				return 0
			neighb1 = np.random.choice(list(neighbs))
			neighbs2 = set(self._adj[neighb1]) & neighbs
			neighbs2 -= {node}
			if not neighbs2:
				return 0
			neighb2 = np.random.choice(list(neighbs2))
		
		co_neighbs = set(self._adj[neighb1]) & set(self._adj[neighb2])
		self.add_e_pair(node, neighb1)
		self.add_e_pair(node, neighb2)
		self.remove_e_pair(neighb1, neighb2)

		if co_neighbs:
			co_neighb = np.random.choice(list(co_neighbs))
			self.add_e_pair(node, co_neighb)

		return 0
	
	def char_squeeout(self, node=None): 
		""" reverse of squeein """
		neighbs = self._adj[node]
		try:
			np.random.shuffle(neighbs)
		except:
			pass
		for neighb1 in neighbs:
			cneighbs = list(self.coneighbours([neighb1, node]))
			if len(cneighbs) > 1:
				self.remove_e_pair(node, neighb1)
				self.add_e_pair(cneighbs[0], cneighbs[1])
				return 0
	
	@classmethod
	def From_Dir(cls, dir, Completeness = 1):
		G = Placoderm.from_edgelist(dir)
		G.graph['completeness'] = Completeness
		G.graph['dir']=dir
		G.Initialise_Symmetry()
		G.attr_from_csv(dir)
		return G

	def attr_from_csv(self, dir):
		""" import character data from csv, add to properties 
		and save property titles """
		attrs = pd.read_csv(dir + "C_Data.txt", sep=" ", index_col=0)
		attr = attrs.to_dict(orient='index')
		nx.set_node_attributes(self, attr)
		return 0

	@classmethod
	def from_edgelist(cls, dir):
		""" derives Morphling from edgelist format """
		file = dir + "G_Data.txt"
		G = nx.read_edgelist(file, nodetype=int)
		G.__class__ = Placoderm
		return G
	

### Business End ###
def main(argv):
	""" Test mutator """	
	G1 = Placoderm.From_Dir(argv[1])
	
	G1.Initialise_Symmetry()
	print(G1.nodes.data())
	G1,a = G1.mutator(move=np.random.choice(range(9)))

	#cProfile.runctx("Profiler(G1)", {'G1':G1, 'Profiler':Profiler}, {})
	G1.draw()
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)