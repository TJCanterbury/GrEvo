#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: all_aligns.sh
# Description: outputs all MI_GRAAL results
# Arguments: none
# Date: 25 April 2021

#### Generate Similarity files for characters/nodes/vertices ####
#Rscript --vanilla char_dists.R ../Data/b.csv ../Data/a.csv ab
#Rscript --vanilla char_dists.R ../Data/a.csv ../Data/c.csv ac
#Rscript --vanilla char_dists.R ../Data/d.csv ../Data/a.csv ad
#Rscript --vanilla char_dists.R ../Data/f.csv ../Data/a.csv af
#Rscript --vanilla char_dists.R ../Data/b.csv ../Data/c.csv bc
#Rscript --vanilla char_dists.R ../Data/d.csv ../Data/b.csv db
#Rscript --vanilla char_dists.R ../Data/d.csv ../Data/c.csv dc
#Rscript --vanilla char_dists.R ../Data/d.csv ../Data/f.csv df
#Rscript --vanilla char_dists.R ../Data/f.csv ../Data/b.csv bf
#Rscript --vanilla char_dists.R ../Data/f.csv ../Data/c.csv cf

#### Run MI-GRAAL for alignments and measuring distance ####
./MI-GRAALRunner.py ../Data/b.gw ../Data/a.gw ../Results/ba -p $1 -q ../Data/ab.txt
./MI-GRAALRunner.py ../Data/a.gw ../Data/c.gw ../Results/ac -p $1 -q ../Data/ac.txt
./MI-GRAALRunner.py ../Data/d.gw ../Data/a.gw ../Results/da -p $1 -q ../Data/ad.txt
./MI-GRAALRunner.py ../Data/f.gw ../Data/a.gw ../Results/fa -p $1 -q ../Data/af.txt
./MI-GRAALRunner.py ../Data/b.gw ../Data/c.gw ../Results/bc -p $1 -q ../Data/bc.txt
./MI-GRAALRunner.py ../Data/d.gw ../Data/b.gw ../Results/db -p $1 -q ../Data/db.txt
./MI-GRAALRunner.py ../Data/d.gw ../Data/c.gw ../Results/dc -p $1 -q ../Data/dc.txt
./MI-GRAALRunner.py ../Data/d.gw ../Data/f.gw ../Results/df -p $1 -q ../Data/df.txt
./MI-GRAALRunner.py ../Data/f.gw ../Data/b.gw ../Results/fb -p $1 -q ../Data/bf.txt
./MI-GRAALRunner.py ../Data/f.gw ../Data/c.gw ../Results/fc -p $1 -q ../Data/cf.txt