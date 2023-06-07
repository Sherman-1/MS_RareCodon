#!/bin/bash


awk 'BEGIN {FS=OFS="\t"} {split($9,a,";"); split(a[1],b,"_"); if(b[5] ~ /^[012]$/) {$8=b[5]}; print $0}' input/mapping_orf_Scer_SGD_noMT.gff > input/test.gff