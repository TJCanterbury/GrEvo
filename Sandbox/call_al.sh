#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: all_aligns.sh
# Description: outputs all MI_GRAAL results
# Arguments: none
# Date: 25 April 2021

#### Generate Similarity files for characters/nodes/vertices ####
Rscript --vanilla char_dists.R 

#### Run MI-GRAAL for alignments and measuring distance ####
./CGRAAL_unix64 ../Data/b.gw ../Data/a.gw ../Data/ab.txt ../Results/ba_numbers.file ../Results/ba_names.file
./CGRAAL_unix64 ../Data/a.gw ../Data/c.gw ../Data/ac.txt ../Results/ac_numbers.file ../Results/ac_names.file
./CGRAAL_unix64 ../Data/d.gw ../Data/a.gw ../Data/ad.txt ../Results/da_numbers.file ../Results/da_names.file
./CGRAAL_unix64 ../Data/f.gw ../Data/a.gw ../Data/af.txt ../Results/fa_numbers.file ../Results/fa_names.file
./CGRAAL_unix64 ../Data/b.gw ../Data/c.gw ../Data/bc.txt ../Results/bc_numbers.file ../Results/bc_names.file
./CGRAAL_unix64 ../Data/d.gw ../Data/b.gw ../Data/db.txt ../Results/db_numbers.file ../Results/db_names.file
./CGRAAL_unix64 ../Data/d.gw ../Data/c.gw ../Data/dc.txt ../Results/dc_numbers.file ../Results/dc_names.file
./CGRAAL_unix64 ../Data/d.gw ../Data/f.gw ../Data/df.txt ../Results/df_numbers.file ../Results/df_names.file
./CGRAAL_unix64 ../Data/f.gw ../Data/b.gw ../Data/bf.txt ../Results/fb_numbers.file ../Results/fb_names.file
./CGRAAL_unix64 ../Data/f.gw ../Data/c.gw ../Data/cf.txt ../Results/fc_numbers.file ../Results/fc_names.file