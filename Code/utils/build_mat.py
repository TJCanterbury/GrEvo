#!/usr/bin/env python3

""" Build parsimony distance matrix from printed incomplete output data """

__appname__ = 'build_mat.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

import sys
from posixpath import sep
import pandas as pd
import os

def find_filenames( path_to_dir, suffix="_final_stats.txt" ):
	""" collect filenames """

	filenames = os.listdir(path_to_dir)
	files = { filename[:-len(suffix)]: (path_to_dir + filename) for filename in filenames if filename.endswith( suffix ) }

	return files

def dist_matrix(G_Data, Pars_Data, file):
	G_files = find_filenames(G_Data, ".txt")
	Data_table = pd.DataFrame(0, columns=G_files, index=G_files)
	Pars_table = pd.read_csv(Pars_Data, sep=" ")
	for index, row in Pars_table.iterrows():
		Data_table.loc[row[0], row[1]] = row[2]
		Data_table.loc[row[1], row[0]] = row[2]
	Data_table.to_csv(file, header=True, index=True)
	return 0

def main(argv):
	dist_matrix(argv[1], argv[2], argv[3])

	return 0


if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)