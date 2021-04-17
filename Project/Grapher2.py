#!/usr/bin/env python3

""" This script can either otuput a plot of a given ID or the AIC results of multiple linear 
and non linear models of temperature performance curves """

__appname__ = 'Mod.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.Godknows'

## imports ##
import sys # module to interface our program with the operating system
import numpy as np
import networkx as nx 
import matplotlib.pyplot as plt

## constants ##
vertices = []
vertices_no = 0
graph = []

## functions ##
def add_vertex(v):
    global graph
    global vertices_no
    global vertices
    if v in vertices:
        print("Vertex ", v, " already exists")
    
    else:
        vertices_no = vertices_no + 1
        vertices.append(v)
        if vertices_no > 1:
            for vertex in graph:
                vertex.append(0)
        temp = []
        for i in range(vertices_no):
            temp.append(0)
        graph.append(temp)

def add_edge(v1, v2, e):
    global graph
    global vertices_no
    global vertices
    if v1 not in vertices:
        print("Vertex ", v1, " doesnt exists")
    elif v2 not in vertices:
        print("Vertex ", v2, " doesnt exists")
    else:
        index1 = vertices.index(v1)
        index2 = vertices.index(v2)
        graph[index1][index2] = e

def print_graph():
    global graph
    global vertices_no

    for i in range(vertices_no):
        for j in range(vertices_no):
            if graph[i][j] != 0:
                print(vertices[i], " -> ", vertices[j], \
                    " edge weight: ", graph[i][j])

def show_graph(adjacency_matrix):
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    gr.add_edges_from(edges)
    nx.draw(gr, node_size=500)
    plt.show()

def main(argv):
    add_vertex("p")
    add_vertex("o1")
    add_vertex("o2")
    add_vertex("g1")
    add_vertex("g2")
    add_vertex("v")
    add_vertex("pi1")
    add_vertex("pi2")
    add_vertex("y1")
    add_vertex("y2")
    add_vertex("c1")
    add_vertex("c2")
    add_vertex("l1")
    add_vertex("l2")

    add_edge('p', 'o1', 1)
    add_edge('p', 'o2', 1)
    add_edge('p', 'v', 1)
    add_edge('g1', 'o1', 1)
    add_edge('g2', 'o2', 1)
    add_edge('v', 'g1', 1)
    add_edge('v', 'g2', 1)
    add_edge('v', 'o1', 1)
    add_edge('v', 'o2', 1)
    add_edge('v', 'c1', 1)
    add_edge('v', 'c2', 1)
    add_edge('v', 'l1', 1)
    add_edge('v', 'l2', 1)
    add_edge('v', 'y1', 1)
    add_edge('v', 'y2', 1)
    add_edge('pi1', 'g1', 1)
    add_edge('pi1', 'c1', 1)
    add_edge('pi1', 'l1', 1)
    add_edge('pi1', 'y1', 1)
    add_edge('pi2', 'g2', 1)
    add_edge('pi2', 'c2', 1)
    add_edge('pi2', 'l2', 1)
    add_edge('pi2', 'y2', 1)
    add_edge('g1', 'c1', 1)
    add_edge('g2', 'c2', 1)
    add_edge('l1', 'c1', 1)
    add_edge('l2', 'c2', 1)
    add_edge('l1', 'y1', 1)
    add_edge('l2', 'y2', 1)


    print_graph()
    print(graph)
    graph1 = np.array(graph)
    show_graph(graph1)
    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)