#!/usr/bin/env python3

""" Build parsimony distance matrix """

__appname__ = 'Pars_GrEv.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
from GrEvo import climber
import pandas as pd
import os
from Morphlings import Plates
import itertools as it

def build_morphling(G_Data, C_Data, Completeness=1):
	""" Builds morphling class with the provided data"""
	G = Plates.from_edgelist(G_Data)
	G.graph['completeness'] = Completeness
	G.attr_from_csv(C_Data)
	return G

def hh_printer(results, best_morphs, best_moves, best_aln):
	""" print pairwise GrEvolutions with the morphs and moves taken and the resulting alignment """
	aln_file = results + "_aln.txt"
	aln = pd.DataFrame.from_dict(best_aln, orient='index')
	aln.to_csv(aln_file)

	Gen= 0
	for i in best_morphs:
			printer = "../Results/" + "morphling_" + results + \
				  "_" + str(Gen) + ".png"
			i.__str__(printer)
			Gen += 1
	
	moves_results = "../Results/" + results + "_moves.txt"
	for i in best_moves:
			with open(moves_results, "a") as myfile:
				myfile.write(i)
	
	return 0

def hh_climb(G1, G2, Result, repeats=10):
	""" Repeat GrEvo algorithm 10 times to get the smallest possible parsomny value"""
	if G1.graph['completeness'] < G2.graph['completeness']:
		G1, G2 = G2, G1
	best_parsimony = 100
	Gen = 0
	while Gen < repeats:
		Gen += 1
		Morphs, moves, aln, parsimony = climber(G1, G2, goal = best_parsimony+1)
		
		if parsimony < best_parsimony:
			best_morphs = Morphs
			best_moves = moves
			best_aln = aln
			best_parsimony = parsimony
	
	hh_printer(Result, best_morphs, best_moves, best_aln)

	return best_parsimony

def find_filenames( path_to_dir, suffix="_final_stats.txt" ):
	""" collect filenames """

	filenames = os.listdir(path_to_dir)
	files = { filename[:-len(suffix)]: (path_to_dir + filename) for filename in filenames if filename.endswith( suffix ) }

	return files

def dist_matrix(G_Data, C_Data, Comp_Data, repeats, file, overwrite=False):
	""" Build parsimony lookup table for use in tree evaluation """
	G_files = find_filenames(G_Data, ".txt")
	Comp_files = find_filenames(Comp_Data, ".txt")
	C_files = find_filenames(C_Data, ".txt")
	Graphs = {}
	
	for G in G_files:
		Graphs[G] = build_morphling(G_files[G], C_files[G])
	
	if os.path.isfile(file):
		Pars_table = pd.read_csv(file)
	else:
		Pars_table = pd.DataFrame(0, columns=Graphs.keys(), index=Graphs.keys())
	print(Pars_table)
	combins = list(it.combinations(G_files.keys(), 2))
	for i in combins:
		if overwrite:
			print(Pars_table)
			print(i)
			score = hh_climb(Graphs[i[0]], Graphs[i[1]], i[0]+i[1], repeats)
			if Pars_table.loc[i[0], i[1]] and score < Pars_table.loc[i[0], i[1]]: # Overwrite old worse scores
				Pars_table.loc[i[0], i[1]] = score
				Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
			
			elif not Pars_table.loc[i[0], i[1]]: # Fill gaps
				Pars_table.loc[i[0], i[1]] = score
				Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
			Pars_table.to_csv(file, header=True)

		elif not Pars_table.loc[i[0], i[1]]:
			print(Pars_table)
			print(i)
			Pars_table.loc[i[0], i[1]] = hh_climb(Graphs[i[0]], Graphs[i[1]], i[0]+i[1], repeats)
			Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
			Pars_table.to_csv(file, header=True)

	return 0

### Business End ###
def main(argv):
	dist_matrix(argv[1], argv[2], argv[3], int(argv[4]), argv[5])
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)