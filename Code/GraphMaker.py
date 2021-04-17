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


## functions ##
def SimVE(size, prob = 0.5):
    ser = np.random.binomial(1, prob, size = size)
    return ser


def AdjacencyMatrix(setV, setE):
	"""generate adjacency matrix for a given morpho character dataset"""
	Adj = np.dot(setV, setE)
	return Adj

def buildgraph(data, cardinality):
    rang = range(cardinality)
    for i in rang: 
        for j in rang: 
            if data[i,j] == 1: 
                G.add_edge(i,j) 



def main(argv):
    cardinality = int(argv[1])
    setV = SimVE(cardinality)
    setE = SimVE(cardinality)

    data = AdjacencyMatrix(setV, setE)
    print(data)
    G = nx.DiGraph()
    buildgraph(data, cardinality)
    nx.draw( G ) 
    plt.show()

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)