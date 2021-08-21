#!/bin/sh
#PBS -lwalltime=12:00:00
#PBS -lselect=1:ncpus=48:mem=64gb

module load anaconda3/personal

python3 ~/Many_GrEvos.py ~/Data/G_Data/ ~/Data/C_Data/ ~/Data/Completeness.csv 1 ~/Steeptest.csv 24 parallel 0 10000

mv Steeptest.csv $HOME
