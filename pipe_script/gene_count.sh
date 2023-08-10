#! /bin/bash

bwa mem -t 4 $1 $2 $3 | samtools view -bS - | samtools sort - > $4
# filter out the reads that are not mapping to the gene set (Flag=4), 
# filter out the reads for the secondary mapping (Flag=256), 
# filter out the chimeric reads (Flag=2048), 
# filter out the reads with an alignment number less than 2
samtools view -F 4 -F 256 -F 2048 $4|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\t'$5'' > $6
