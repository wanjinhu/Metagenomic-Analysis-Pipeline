#! /bin/bash

prodigal -p meta -a $1 -m -d $2 -o $3 -f gff -s $4 -i $5

## demo
# prodigal基因预测
# prodigal -p meta -a DBYX1FB_prot.faa -m -d DBYX1FB_nucl.fna -o DBYX1FB_genes.gff -f gff -s DBYX1FB.stat -i DBYX1FB.contigs_500.fa