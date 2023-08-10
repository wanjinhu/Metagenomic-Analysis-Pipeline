#! /bin/bash

zcat $1 $2|metaphlan --input_type fastq --bowtie2out $3 --output_file $4 --nproc 8
