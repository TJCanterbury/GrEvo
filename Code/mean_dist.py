#!/usr/bin/env python3

""" Finds mean edge correctness across the distance matrices of the different batches of MI-GRAAL results """

__appname__ = 'mean_dist.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import numpy as np
from numpy import genfromtxt
from Dists import find_filenames

## functions ##
def mean_np(files, n, path = "../Results/"):
    """ Finds mean edge correctness across the distance matrices of the files """
    
    data = np.zeros((5, 5))

    for i in range(n):
        file = path + files[i]
        print(file)
        data += genfromtxt(file, delimiter=',')
        print(data)
    data /= n

    return data

def main(argv):
    """ Feeds all fistance matrix files into mean_np and saves the resulting matrix """
    files = find_filenames(path_to_dir = "../Results/", suffix=".csv")
    n = len(files)
    Data = mean_np(files, n)

    np.savetxt("../Results/mean_dist.csv", Data, delimiter=',')

    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)