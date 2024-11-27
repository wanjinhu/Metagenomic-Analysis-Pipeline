#! /bin/bash

zcat $1 $2|metaphlan --input_type fastq --bowtie2out $3 --output_file $4 --nproc 16 --ignore_usgbs --index mpa_vJun23_CHOCOPhlAnSGB_202307

## demo
# 基于metaphlan物种注释
# zcat DBYX1FB.unmap.1.fastq.gz DBYX1FB.unmap.2.fastq.gz|metaphlan --input_type fastq --bowtie2out DBYX1FB_bowtie2.bz2 --output_file DBYX1FB_metaphlan.tsv --nproc 8
# 多个样本合并物种信息
# merge_metaphlan_tables.py *.txt > merged_abundance_table.txt