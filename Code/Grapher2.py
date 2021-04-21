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
from numpy import genfromtxt

## constants ##
vertices = []
vertices_no = 0
graph = []

## functions ##
def unique_labels(np_array_data):
    return np.unique(np_array_data)

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

def add_all_v(data):
    Uniques = unique_labels(data)
    for v in Uniques:
        add_vertex(Uniques[v])

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

def add_all_e(data):
    add_edge()

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
    my_data = genfromtxt(argv[1], delimiter=' ')
    #add_all_v(my_data) 
    
    print_graph()
    print(graph)
    graph1 = np.array(graph)
    show_graph(graph1)
    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)