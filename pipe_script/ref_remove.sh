#! /bin/bash

bowtie2 -x $1 -1 $2 -2 $3 -S $4 -p 32  2>$5
samtools fastq -@ 32 -f 4 $4 -1 $6 -2 $7 -s $8

## demo
# 去宿主
# bowtie2 -x /root/database/kneaddata_human/hg37dec_v0.1 
#         -1 DBYX1FB_clean.1.fastq.gz 
#         -2 DBYX1FB_clean.2.fastq.gz 
#         -S DBYX1FB.sam  2>DBYX1FB.mapping.log
# 提取未比对到宿主的序列, -f 4表示：过滤掉未比对到参考基因组的reads（Flag=4）
# samtools fastq -@ 8 -f 4 DBYX1FB.sam -1 DBYX1FB.unmap.1.fastq.gz -2 DBYX1FB.unmap.2.fastq.gz -s DBYX1FB.unmap.single.fastq.gz

## demo
# 去两个宿主序列
# 去宿主1
# bowtie2 -x /root/database/hg38_GCF_000001405.40/GCF_000001405.40/hg38 \
#         -1 B-NJBIAN-10_1.fq.gz \
#         -2 B-NJBIAN-10_2.fq.gz \
#         -S process.sam \
#         -p 32
# 提取未比对到宿主的序列, -f 4表示：过滤掉未比对到参考基因组的reads（Flag=4）
# samtools fastq -@ 32 -f 4 process.sam -1 process_1.fastq.gz -2 process_2.fastq.gz -s process_single.fastq.gz
# 去宿主2
# bowtie2 -x /root/database/Canis_GCF_000002285.5/Canis_GCF_000002285_5 \
#         -1 process_1.fastq.gz \
#         -2 process_2.fastq.gz \
#         -U process_single.fastq.gz \
#         -S B-NJBIAN-10.sam \
#         -p 32
# samtools fastq -@ 32 -f 4 B-NJBIAN-10.sam -1 B-NJBIAN-10.unmap.1.fastq.gz -2 B-NJBIAN-10.unmap.2.fastq.gz -s B-NJBIAN-10.unmap.single.fastq.gz

