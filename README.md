# Welcome to Grevo (graph evolution) repository

### Disclaimer 
MI-GRAAL, MI-GRAALRunner.py, list2leda, ncount, portrait_divergence.py, magnapp_cli_linux64 and much of the UPGMA function are not my own but can be found online readily.

## Directory Structure:
 
 Directory       | Description
 ------------- | -------------
 mi-graal_UPGMA | MI-GRAAL is run with 'bash all_aligns.sh ../pond/' to generate a distance matric and from this a UPGMA evaluated phylogeny is generated.
 magna_UPGMA | magna is run with 'bash all_aligns.sh ../pond/' to generate a distance matric and from this a UPGMA evaluated phylogeny is generated.
 GrEvo | Used to randomize and evolve graphs from one topology to another. Will be used to measure parsimony.

## Requirements:
### System:
Tested on Ubuntu 20.04.2 LTS

### Languages:
- Python 3.8.5
    Imports:
    - Networkx
    - NumPy
    - Pandas
    - itertools
    - os
    - Biopython
    - matplotlib
    - io
    - re
    - sys
- R 3.6.3
    Imports:
    - tidyverse
