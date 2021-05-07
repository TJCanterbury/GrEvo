# Welcome to Grevo (graph evolution) repository

## Directory Structure:

Sandbox is where all my junk code went and in the Data and pond directory you will find many graphs.  

### In Code you shall find:
 
 File       | Description
 ------------- | -------------
 all_aligns.sh | MI-GRAAL is run with 'bash all_aligns.sh ../pond/' to generate a distance matric and from this a UPGMA evaluated phylogeny.
 grevo.py | This python script is used with graphs (eg a.txt) in Data to randomize and evolve graphs. This shall be adapted with a hill climbing algorithm to evolve one graph into another, find the shortest path for this evolution and use the length of this path as a measure of parsimony in future tree evaluations. example usage: 'python3 grevo.py ../Data/a.txt ../Data/b.txt 10'