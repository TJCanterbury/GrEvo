#!/usr/bin/env python3

""" A Quick Implementation of UPGMA (Unweighted Pair Group Method with Arithmetic Mean) """

__appname__ = 'UPGMA.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import re
import os
import numpy as np
import pandas as pd
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt
import networkx as nx 

## Functions ##
# lowest_cell:
#   Locates the smallest cell in the table
def lowest_cell(table):
    # Set default to infinity
    min_cell = float("inf")
    x, y = -1, -1

    # Go through every cell, looking for the lowest
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min_cell:
                min_cell = table[i][j]
                x, y = i, j

    # Return the x, y co-ordinate of cell
    return x, y


# join_labels:
#   Combines two labels in a list of labels
def join_labels(labels, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # Join the labels in the first index
    labels[a] = "(" + labels[a] + "," + labels[b] + ")"

    # Remove the (now redundant) label in the second index
    del labels[b]


# join_table:
#   Joins the entries of a table on the cell (a, b) by averaging their data entries
def join_table(table, a, b):
    # Swap if the indices are not ordered
    if b < a:
        a, b = b, a

    # For the lower index, reconstruct the entire row (A, i), where i < A
    row = []
    for i in range(0, a):
        row.append((table[a][i] + table[b][i])/2)
    table[a] = row
    
    # Then, reconstruct the entire column (i, A), where i > A
    #   Note: Since the matrix is lower triangular, row b only contains values for indices < b
    for i in range(a+1, b):
        table[i][a] = (table[i][a]+table[b][i])/2
        
    #   We get the rest of the values from row i
    for i in range(b+1, len(table)):
        table[i][a] = (table[i][a]+table[i][b])/2
        # Remove the (now redundant) second index column entry
        del table[i][b]

    # Remove the (now redundant) second index row
    del table[b]


# UPGMA:
#   Runs the UPGMA algorithm on a labelled table
def UPGMA(table, labels):
    # Until all labels have been joined...
    while len(labels) > 1:
        # Locate lowest cell in the table
        x, y = lowest_cell(table)

        # Join the table on the cell co-ordinates
        join_table(table, x, y)

        # Update the labels accordingly
        join_labels(labels, x, y)

    # Return the final label
    return labels[0]

def main(argv): 
    # Test table data and corresponding labels
    M_labels = ['Wuttagoonaspis', 'Romundina', 'Brindabellaspis', 'Eurycaraspis', 'Entelognathus']  
    print(M_labels) #A through G
    M = np.loadtxt(open(argv[1], "rb"), delimiter=",")
    M = np.tril(M).tolist()
    
    for j in range(0,5):
        for i in range(0,5):
            M[i] = list(filter(lambda a: a !=0, M[i]))
    print(M)

    tree = UPGMA(M, M_labels)

    handle = StringIO(tree)
    tree = Phylo.read(handle, "newick")
    print(tree)
    tree.ladderize()
    Phylo.draw(tree)

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)