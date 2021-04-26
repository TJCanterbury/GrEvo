#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: Laign.sh
# Description: outputs all L_GRAAL results, issue: with segmentation fault
# Arguments: none
# Date: 26 April 2021

#### Generate degree signatures and character similarities
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
#
#./ncount4.exe ../Data/a.gw ../Data/a
#./ncount4.exe ../Data/b.gw ../Data/b
#./ncount4.exe ../Data/c.gw ../Data/c
#./ncount4.exe ../Data/d.gw ../Data/d
#./ncount4.exe ../Data/f.gw ../Data/f

./L-GRAAL.exe -Q ../Data/b.gw -T ../Data/a.gw -q ../Data/b.ndump2 -t ../Data/a.ndump2 -o ../Results/ba -B ../Data/ab.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/a.gw -T ../Data/c.gw -q ../Data/a.ndump2 -t ../Data/c.ndump2 -o ../Results/ac -B ../Data/ac.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/d.gw -T ../Data/a.gw -q ../Data/d.ndump2 -t ../Data/a.ndump2 -o ../Results/da -B ../Data/ad.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/f.gw -T ../Data/a.gw -q ../Data/f.ndump2 -t ../Data/a.ndump2 -o ../Results/fa -B ../Data/af.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/b.gw -T ../Data/c.gw -q ../Data/b.ndump2 -t ../Data/c.ndump2 -o ../Results/bc -B ../Data/bc.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/d.gw -T ../Data/b.gw -q ../Data/d.ndump2 -t ../Data/b.ndump2 -o ../Results/db -B ../Data/db.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/d.gw -T ../Data/c.gw -q ../Data/d.ndump2 -t ../Data/c.ndump2 -o ../Results/dc -B ../Data/dc.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/d.gw -T ../Data/f.gw -q ../Data/d.ndump2 -t ../Data/f.ndump2 -o ../Results/df -B ../Data/df.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/f.gw -T ../Data/b.gw -q ../Data/f.ndump2 -t ../Data/b.ndump2 -o ../Results/fb -B ../Data/bf.txt -I 10000 -L 36000
./L-GRAAL.exe -Q ../Data/f.gw -T ../Data/c.gw -q ../Data/f.ndump2 -t ../Data/c.ndump2 -o ../Results/fc -B ../Data/cf.txt -I 10000 -L 36000