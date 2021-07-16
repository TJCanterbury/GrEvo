import pandas as pd
import itertools as it
Pars_table1 = pd.read_csv("pars3.csv")
Pars_table2 = pd.read_csv("pars4.csv")

combins = list(it.combinations(list(Pars_table1.columns), 2))

for i in combins:
	
	print(((Pars_table1.loc[i[0], i[1]] - Pars_table2.loc[i[0], i[1]]) / Pars_table1.loc[i[0], i[1]]) * 100)
	if Pars_table1.loc[i[0], i[1]] > Pars_table2.loc[i[0], i[1]]:
		Pars_table1.loc[i[0], i[1]] = Pars_table2.loc[i[0], i[1]]
		Pars_table1.loc[i[1], i[0]] = Pars_table2.loc[i[0], i[1]]
		

print(Pars_table1)

Pars_table1.to_csv("pars7.csv", header=True, index=True)