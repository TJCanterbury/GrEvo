# GrEvo the Graph Evolver is a work in progress

## Use:

Enter the 2 graphs, the one you want to start at and the one you want to evolve to, how many perturbations of the first graph to try each generation of the climbing algorithm, a results file name for storing the steps taken... if it ever completes and the number of attempts before a graph is determined a dead end and finally the max number of steps you'll allow (for instance if you know it can be done in 4 moves you could set the limit to 4 moves).

to run all graph comparisons run the following from the Code directory:
python3 Pars_GrEv.py ../Data/G_Data/ ../Data/C_Data/ ../Data/Completeness.csv 1 ../Results/test.csv 10 parallel 0