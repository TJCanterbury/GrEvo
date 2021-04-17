#!/bin/bash
# Author: Tristan JC tjc19@ic.ac.uk
# Script: CompileLaTeX.sh
# Description: compiles tex file giving you a pdf
# Arguments: 1 tex file
# Date: 8 Oct 2020
if [ $# -ne 1 ]
    then
        echo "missing or too many arguments, try again with a single .tex file"
        exit
fi
if [ ${1: -4} != ".tex" ]
    then
        read -p "$1 is not labelled as a .tex file, try again with a single .tex file"
        exit
fi

x=$1
y=${x%.tex}
aspell -t -c $y.tex
pdflatex $y.tex
pdflatex $y.tex
bibtex $y
pdflatex $y.tex
pdflatex $y.tex
mv $y.pdf ../Results
evince ../Results/$y.pdf &

## Cleanup

rm *~
rm *.aux
rm *.dvi
rm *.log
rm *.nav
rm *.out
rm *.snm
rm *.toc
rm *.bbl
rm *.blg
rm *.fls
rm *.fdb*
rm *.Rout
rm *concordance*
rm .RData