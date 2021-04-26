#!/usr/bin/env python3

""" Uses generic networkx code to generate graphs from our data """

__appname__ = 'nxcode.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx 
import matplotlib.pyplot as plt
import csv

## functions ##
def readx(filex):
    with open(filex) as infile:
        csv_reader = csv.reader(infile, delimiter=' ')
        G = nx.Graph(csv_reader)
    return G

def drawx(G):
    nx.draw_networkx(G, with_labels = True)
    plt.show()

def main(argv):
    G1 = readx(argv[1])
    drawx(G1)

    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)