#!/usr/bin/env python3

""" Merge matrices to find smallest scores """

__appname__ = 'merge.py'
__author__ = 'Tristan JC (tjc19@ic.ac.uk)'
__version__ = '0.0.1'

## imports ##
import sys
import pandas as pd
import itertools as it

def merge(M1, M2, name):
	Pars_table1 = pd.read_csv(M1)
	Pars_table2 = pd.read_csv(M2)

	combins = list(it.combinations(list(Pars_table1.columns), 2))

	for i in combins:
		
		print(((Pars_table1.loc[i[0], i[1]] - Pars_table2.loc[i[0], i[1]]) / Pars_table1.loc[i[0], i[1]]) * 100)
		if Pars_table1.loc[i[0], i[1]] == 0 and Pars_table2.loc[i[0], i[1]] != 0:
			Pars_table1.loc[i[0], i[1]] = Pars_table2.loc[i[0], i[1]]
			Pars_table1.loc[i[1], i[0]] = Pars_table2.loc[i[0], i[1]]

		if Pars_table1.loc[i[0], i[1]] > Pars_table2.loc[i[0], i[1]] and Pars_table2.loc[i[0], i[1]] != 0:
			Pars_table1.loc[i[0], i[1]] = Pars_table2.loc[i[0], i[1]]
			Pars_table1.loc[i[1], i[0]] = Pars_table2.loc[i[0], i[1]]
		
			

	print(Pars_table1)

	Pars_table1.to_csv(name, header=True, index=True)

	return 0


### Business End ###
def main(argv):
	merge(argv[1], argv[2], argv[3])

	return 0

if __name__ == "__main__": 
	"""Makes sure the "main" function is called from command line"""  
	status = main(sys.argv)
	sys.exit(status)