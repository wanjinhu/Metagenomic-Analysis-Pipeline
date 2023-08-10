#! /bin/bash

bowtie2 -x $1 -1 $2 -2 $3 -S $4  2>$5
# -f 4 means retain the reads that are not mapping to the reference genome 
samtools fastq -@ 8 -f 4 $4 -1 $6 -2 $7 -s $8
