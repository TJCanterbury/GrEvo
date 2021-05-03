#!/usr/bin/env python3

""" Uses generic networkx code to generate graphs from our data """

__appname__ = 'mean_dist.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import re
import os
import numpy as np
from numpy import genfromtxt
import pandas as pd
import networkx

## functions ##
def find_filenames( path_to_dir = "../Results/", suffix=".csv" ):
    filenames = os.listdir(path_to_dir)
    files = [ filename for filename in filenames if filename.endswith( suffix ) ]
    print(files)
    return files

def mean_np(files, n, path = "../Results/"):
    data = np.zeros((5, 5))

    for i in range(n):
        file = path + files[i]
        print(file)
        data += genfromtxt(file, delimiter=',')
        print(data)
    data /= n

    return data

def main(argv):

    files = find_filenames()
    n = len(files)
    Data = mean_np(files, n)

    np.savetxt("../Results/mean_dist.csv", Data, delimiter=',')

    return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)