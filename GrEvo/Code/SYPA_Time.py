#!/usr/bin/env python3

""" Generate SYPA alignment scores as distance matrix for analysis of performance """

__appname__ = 'SYPA_Scores.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
from SYPA import SYPA
import sys
import pandas as pd
import os
from Morphlings import Plates
import itertools as it
from joblib import Parallel, delayed
from Morphlings import Plates
import numpy as np
import matplotlib.pyplot as plt

from numpy.polynomial.polynomial import polyfit

### Functions ###
def build_morphling(G_Data, C_Data, Completeness):
	""" Builds morphling class with the provided data"""
	G = Plates.from_edgelist(G_Data)
	G.graph['completeness'] = Completeness
	G.attr_from_csv(C_Data)
	return G

def find_filenames( path_to_dir, suffix="_final_stats.txt" ):
	""" collect filenames """

	filenames = os.listdir(path_to_dir)
	files = { filename[:-len(suffix)]: (path_to_dir + filename) for filename in filenames if filename.endswith( suffix ) }

	return files

def SYPA_job(G1, G2, repeats=1000):
	""" Repeat GrEvo algorithm 10 (for example) times to get the smallest possible parsimony value"""
	if G1.graph['completeness'] < G2.graph['completeness'] or G1.number_of_nodes() > G2.number_of_nodes():
		G1, G2 = G2, G1
	
	CS3, aln = SYPA(G1, G2, repeat=repeats)
	print(CS3)
	return CS3

def dist_matrix(G_Data, C_Data, Comp_Data, file, num_cores):
	""" Build parsimony lookup table for use in tree evaluation """
	G_files = find_filenames(G_Data, ".txt")
	Completeness = pd.read_csv(Comp_Data,index_col=0,header=None)
	C_files = find_filenames(C_Data, ".txt")
	Graphs = {}
	
	for G in G_files:
		Graphs[G] = build_morphling(G_files[G], C_files[G], int(Completeness.loc[G]))
	
	combins = list(it.combinations(G_files.keys(), 2))

	results = np.empty((len(combins),100))

	for j in range(100):
		results[:,j] = Parallel(n_jobs=num_cores)(delayed(SYPA_job)(G1=Graphs[i[0]], G2=Graphs[i[1]],\
			repeats=j+1) for i in combins)
	sum_of_rows = results.sum(axis=1)

	results = results / sum_of_rows[:, np.newaxis]
	results= pd.DataFrame(results).reset_index().melt("index")
	print(results)
	np.savetxt("GenvS3_SYPA.txt", results)
	plt.scatter(results["variable"]+1, results["value"])
	
	plt.show()
	return 0

### Business End ###
def main(argv):
	dist_matrix(G_Data=argv[1], C_Data=argv[2], Comp_Data=argv[3], \
		file=argv[4], num_cores=int(argv[5]))
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)