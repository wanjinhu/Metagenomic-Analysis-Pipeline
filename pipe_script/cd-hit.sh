#! /bin/bash

cd-hit -i $1 -o $2 -c 0.95 -T 8 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
grep '>' $2|awk -F ' ' '{print $1}'|sed 's/>//g' > $3
seqtk subseq $4 $3 > $5
bwa index $5 -p $6
echo -e "gene\tlength" > $7
bioawk -c fastx '{print $name, length($seq)}' $5 >> $7
