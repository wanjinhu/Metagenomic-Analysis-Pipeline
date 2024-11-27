#! /bin/bash

# /root/biosoft/cd-hit-v4.8.1-2019-0228/cd-hit
cd-hit -i $1 -o $2 -c 0.95 -T 0 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
grep '>' $2|awk -F ' ' '{print $1}'|sed 's/>//g' > $3
/root/miniconda3/bin/seqtk subseq $4 $3 > $5

# /root/biosoft/bwa/bwa
bwa index $5 -p $6
echo -e "gene\tlength" > $7
# /root/biosoft/bioawk/bioawk
bioawk -c fastx '{print $name, length($seq)}' $5 >> $7

## demo
## cd-hit去冗余,先用蛋白序列去冗余
# 由于不同的核酸序列翻译后可能产生相同的蛋白质序列，因此对核酸序列的去冗余需要基于蛋白质序列
# 基于蛋白质序列文件与核酸序列文件中序列号的对应关系，根据去冗余后的蛋白质序列文件的序列号筛选核酸序列
# cd-hit -i prot.faa -o prot_nonerude.faa -c 0.95 -T 8 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
# less prot_nonerude.faa|grep '>'|awk -F ' ' '{print $1}'|sed 's/>//g' > prot_nonerude.list
# seqtk subseq nucl.fna prot_nonerude.list > nucl_nonerude.fna
## 非冗余基因集构建bwa索引,后续样本去比对
# bwa index nucl_nonerude.fna -p geneset_bwa
# 计算非冗余基因集中各个基因的长度，为了后面可能计算基因的RPKM等信息，需要基因的长度信息
# (echo -e "gene\tlength"; bioawk -c fastx '{print $name, length($seq)}' nucl_nonerude.fna) > geneset_length.txt