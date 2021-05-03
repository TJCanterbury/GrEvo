#!/usr/bin/env python3

""" This script will be the start of my attempt to simulate phenotypic evolution through a hill climbing algorithm applied to perturbations of graphs.
These vertices of the graph represent morphological characters and edges represent their physical connections with each other. The hypothesis is 
that as connection are based on size and position of characters, how these edges change will be an effective model of phenotypic evolution.

To test this model I will use this code to find the least number of changes needed to go from one anatomical network to the next, using edge correctness -- 
estimated with MI-GRAAL -- to measure distance in isomorphism between the anatomical networks. These changes will then be the most parsimonious explanations for 
how one species may be translated into another and so from there we can build a tree, where the most parsimonious translations are 
Then I will build a phylogeny based on these events. """

__appname__ = 'grevo.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys # module to interface our program with the operating system
import networkx as nx 
import matplotlib.pyplot as plt
import csv

