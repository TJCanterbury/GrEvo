#!/usr/bin/env python3

""" Build parsimony distance matrix """

__appname__ = 'Many_GrEvos.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
from Steepest_GrEvo import Steep_GrEv
import sys
import pandas as pd
import os
from Morphlings import Plates
import itertools as it
from joblib import Parallel, delayed
from Morphlings import Plates

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

def GrEvo_job(G1, G2, G1_name, G2_name, Climb_repeats, goal = 1000, breadth=10000):
	""" Repeat GrEvo algorithm 10 (for example) times to get the smallest possible parsimony value"""
	if G1.graph['completeness'] < G2.graph['completeness'] or G1.number_of_nodes() > G2.number_of_nodes():
		G1, G2 = G2, G1
	if goal == 0:
		goal = 1000
	best_parsimony = goal
	Gen = 0
	while Gen < Climb_repeats:
		Gen += 1
		temp, parsimony = Steep_GrEv(G1, G2, G1_name=G1_name, G2_name=G2_name, goal = best_parsimony, Breadth=breadth)
		
		print(Gen)
		if parsimony < best_parsimony:
			best_parsimony = parsimony
			
	# if both graphs are assumed complete repeat the process but evolving in the reverse direction
	if G1.graph['completeness'] + G2.graph['completeness'] == 2:
		G1, G2 = G2, G1
		
		Gen = 0
		while Gen < Climb_repeats:
			Gen += 1
			temp, parsimony = Steep_GrEv(G1, G2, G1_name=G1_name, G2_name=G2_name, goal = best_parsimony, Breadth=breadth)

			if parsimony < best_parsimony:
				best_parsimony = parsimony
	
	print(G1_name + " X " + G2_name + ": " + str(best_parsimony))
	return [G1_name, G2_name, best_parsimony]

def dist_matrix(G_Data, C_Data, Comp_Data, repeats, file, num_cores, mode="simple", overwrite=False, breadth=10000):
	""" Build parsimony lookup table for use in tree evaluation """
	G_files = find_filenames(G_Data, ".txt")
	Completeness = pd.read_csv(Comp_Data,index_col=0,header=None)
	C_files = find_filenames(C_Data, ".txt")
	Graphs = {}
	
	for G in G_files:
		Graphs[G] = build_morphling(G_files[G], C_files[G], int(Completeness.loc[G]))
	
	if os.path.isfile(file):
		Pars_table = pd.read_csv(file)
	else:
		Pars_table = pd.DataFrame(0, columns=Graphs.keys(), index=Graphs.keys())

	combins = list(it.combinations(G_files.keys(), 2))

	if mode == "parallel":
		print(Pars_table)
		
		if overwrite:
			Results = Parallel(n_jobs=num_cores)(delayed(GrEvo_job)(G1=Graphs[i[0]], G2=Graphs[i[1]], G1_name = i[0], \
				G2_name=i[1], Climb_repeats=repeats, goal=Pars_table.loc[i[0], i[1]], breadth=breadth) for i in combins)
			for i in Results:
				if Pars_table.loc[i[0], i[1]] > i[2]:
					Pars_table.loc[i[0], i[1]] = i[2]
					Pars_table.loc[i[1], i[0]] = i[2]
		else:
			Results = Parallel(n_jobs=num_cores)(delayed(GrEvo_job)(G1=Graphs[i[0]], G2=Graphs[i[1]], G1_name = i[0], \
				G2_name=i[1], Climb_repeats=repeats, breadth=breadth) for i in combins if not Pars_table.loc[i[0], i[1]])
			for i in Results:
				Pars_table.loc[i[0], i[1]] = i[2]
				Pars_table.loc[i[1], i[0]] = i[2]
		
		Pars_table.to_csv(file, header=True, index=True)
		print(Pars_table)
		os.system('play -nq -t alsa synth {} sine {}'.format(0.1, 150))

	if mode == "simple":
		for i in combins:
			if overwrite:
				print(Pars_table)
				print(i)
				score = GrEvo_job(Graphs[i[0]], Graphs[i[1]], i[0], i[1], repeats)[2]
				if Pars_table.loc[i[0], i[1]] and score < Pars_table.loc[i[0], i[1]]: # Overwrite old worse scores
					Pars_table.loc[i[0], i[1]] = score
					Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
				
				elif not Pars_table.loc[i[0], i[1]]: # Fill gaps
					Pars_table.loc[i[0], i[1]] = score
					Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
				Pars_table.to_csv(file, header=True, index=True)

			elif not Pars_table.loc[i[0], i[1]]:
				print(Pars_table)
				print(i)
				result = GrEvo_job(Graphs[i[0]], Graphs[i[1]], i[0], i[1], repeats)
				Pars_table.loc[i[0], i[1]] = result[2]
				Pars_table.loc[i[1], i[0]] = Pars_table.loc[i[0], i[1]]
				Pars_table.to_csv(file, header=True, index=True)

	return 0

### Business End ###
def main(argv):
	dist_matrix(G_Data=argv[1], C_Data=argv[2], Comp_Data=argv[3], \
		repeats=int(argv[4]), file=argv[5], num_cores=int(argv[6]), mode=argv[7], overwrite=int(argv[8]),
		breadth=int(argv[9]))
	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)