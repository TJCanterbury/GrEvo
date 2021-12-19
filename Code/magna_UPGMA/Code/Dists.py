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
        labels=re.findall('[A-Z][^A-Z]*', label)
        l1 = labels[0]
        l2=labels[1].split("_")[0]
        filep = "../Results/" + file
        s = get_s(filep)
        data.loc[i-1, 0] = l1
        data.loc[i-1, 1] = l2
        data.loc[i-1, 2] = s
        i += 1
    print(data)
    return data

def make_mat(data):
    """ Turn list into adjacency matrix for use in UPGMA """
    
    data = data.pivot(index=0,columns=1)

    return data

def main(argv):
    """ Produce distance matrix from the MAGNA derived edge correctness """
    
    path = argv[1] # path to outputs of MAGNA
    files = find_filenames(path)
    data = list_vals(files)
    data.to_csv(argv[2])
    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)