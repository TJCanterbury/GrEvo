#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: MI_Runner.sh
# Description: outputs MI_GRAAL results for a given network pair
# Arguments: 2 network files with .gw leda format
# Date: 21 April 2021

### Check arguments ###
if [ $# -ne 4 ]
    then
        echo "missing or too many arguments, try again with a 2 .gw network files"
        exit
fi

./MI-GRAALRunner.py $1 $2 $3 -p $4