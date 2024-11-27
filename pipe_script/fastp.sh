#! /bin/bash

/root/miniconda3/bin/fastp -i $1 -o $2 -I $3 -O $4 -w 32 -h $5 -j $6

## demo
# fastp -i DBYX1FB_1.fastq.gz \
#       -o DBYX1FB_clean.1.fastq.gz \
#       -I DBYX1FB_2.fastq.gz \
#       -O DBYX1FB_clean.2.fastq.gz \
#       -w 8 -h DBYX1FB.html -j DBYX1FB.json
