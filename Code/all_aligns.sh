#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: all_aligns.sh
# Description: outputs all MI_GRAAL results
# Arguments: none
# Date: 25 April 2021

#### Generate Similarity files for characters/nodes/vertices ####
Rscript --vanilla char_dists.R 16

rm ../Results/mean_dist.csv

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
    #### Run MI-GRAAL for alignments and measuring distance ####
    sh GRAAL_it.sh

    python3 /home/tristan/Documents/Masters/Code/Dists.py ../Results/Dist_matrix_$i.csv

done
cat GRAAL_it.sh
rm GRAAL_it.sh
python3 mean_dist.py
python3 UPGMA.py ../Results/mean_dist.csv 