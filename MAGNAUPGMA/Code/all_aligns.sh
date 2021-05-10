#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: all_aligns.sh
# Description: outputs all MAGNA results
# Arguments: none
# Date: 25 April 2021

echo Hello, Where are the graphs and similarity scores located?
read path

echo "Would you like to specify th parameters? (y/n)"
read specify

if [ $specify = y ] ; then
    echo From a scale of 0 to 1 how accurate are those node similarities?
    read alpha
    echo How many alighments should we align per generation?
    read population
    echo How many Generations?
    read generations
    echo Over how many threads should the process be parralized?
    read threads
fi

if [ $specify = n ] ; then
    alpha="0.5"
    population="1000"
    generations="1000"
    threads="10" 
fi

#### Generate Similarity files for characters/nodes/vertices ####
Rscript --vanilla char_dists.R $path $alpha $population $generations $threads

#### Build alignments and output scores #### 
sh MAGNA_it.sh

#### Build Distance matrix for taxa from alignment, WARNING: If more taxa are added this will need editing!!
python3 /home/tristan/Documents/Masters/Code/Dists.py ../Results/ ../Results/Dist_matrix.csv

#### Build tree from distance matrix using UPGMA, save the tree yourself ####
python3 UPGMA.py ../Results/Dist_matrix.csv