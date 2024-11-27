#! /bin/bash

mkdir $4
salmon quant -i $1 -l A -p 16 --meta -1 $2 -2 $3 -o $4/$5
# 获取salmon的TPM结果
less $4/$5/quant.sf|cut -f1,4|sed '1d'|sed '1i gene\t'$6'' > $7
# 获取salmon的reads counts结果
less $4/$5/quant.sf|cut -f1,5|sed '1d'|sed '1i gene\t'$6'' > $8

# demo
# salmon index -t nucl_nonerude.fna -p 16 -i salmon_index
# mkdir out_salmon
# salmon quant -i salmon_index -l A -p 16 --meta -1 test.unmap.1.fastq.gz -2 test.unmap.2.fastq.gz -o out_salmon/test.quant
# 获取salmon的TPM结果
# less out_salmon/test.quant/quant.sf|cut -f1,4|sed '1d'|sed '1i gene\tbowtie2' > salmon_TPM.count
# 获取salmon的reads counts结果
# less out_salmon/test.quant/quant.sf|cut -f1,5|sed '1d'|sed '1i gene\tbowtie2' > salmon_reads.count
