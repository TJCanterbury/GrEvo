#!/usr/bin/env python3

""" Produce distance matrix from the MI-GRAAL edge correctness output from .results files, including all combinations """

__appname__ = 'Dists.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import re
import os
import numpy as np
import pandas as pd
import networkx

## functions ##
def get_s(file):
    """ retrieve edge correctness from file of a given name/path """
    s = pd.read_csv(file, delimiter = " ", index_col=0) 

    s = float(s.iloc[4]) 
    print(s)
    return(s)

def find_filenames( path_to_dir, suffix="_final_stats.txt" ):
    """ collect filenames """

    filenames = os.listdir(path_to_dir)
    files = [ filename for filename in filenames if filename.endswith( suffix ) ]
    print(files)
    return files

def list_vals(files):
    """ make data structure for distances and assign them """
    data = np.ones((len(files), 3))
    data = pd.DataFrame(data)
    i = 1
    for file in files:
        label = os.path.splitext(file)[0]
        filep = "../Results/" + file
        s = get_s(filep)
        data.loc[i-1, 0] = label[0]
        data.loc[i-1, 1] = label[1]
        data.loc[i-1, 2] = 1 - s
        i += 1
    print(data)
    return data

def make_mat(data):
    """ Turn list into adjacency matrix for use in UPGMA """

    letters = []
    for i in range(6):
        letter = ['a', 'b', 'c', 'd', 'e', 'f'][i]
        letters.append(letter)
    G = networkx.OrderedGraph()
    edgeList = data.values.tolist()

    G.add_nodes_from(letters)
    for i in range(len(edgeList)):
        G.add_edge(edgeList[i][0], edgeList[i][1], weight=edgeList[i][2])
    A = networkx.adjacency_matrix(G).A
    print(G.nodes)
    return A

def main(argv):
    """ Produce distance matrix from the MAGNA derived edge correctness """
    
    path = argv[1] # path to outputs of MAGNA
    files = find_filenames(path)
    data = list_vals(files)
    matrix = make_mat(data)
    print(matrix)
    np.savetxt(argv[2], matrix, delimiter=',') # argv[2] is the file location you want to put the distance matrix

    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)