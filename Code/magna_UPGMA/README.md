# Using MAGNA to align anatomical networks to measure phylogenetic distance between placoderm roof structures, UPGMA used to build Phylogeny

## Use:

Enter the code directory and run the following:
//start//
bash all_aligns.sh 
//end//

You will then be asked to provide a path to where the graphs and their node character data can be found, enter ../Data/
Then it will ask if you want to specifiy the arguments yourself, if you enter n default arguments will be used (WARNING: default of 10 threads, check you have this cpu capacity, if you have less than 10, answer y and specify your own thread count.)