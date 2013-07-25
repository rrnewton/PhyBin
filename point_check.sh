#!/bin/bash

arg1=$1
arg2=$2

phybin -p2 -b0.02 -n10 --complete --editdist=0 --rfdist Wolbachia/trees/*.$arg1.codon Wolbachia/trees/*.$arg2.codon
