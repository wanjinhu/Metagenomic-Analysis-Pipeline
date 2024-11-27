#! /bin/bash

bowtie2 -x $1 -1 $2 -2 $3 -S $4 -p 16
samtools view $4 -bS | samtools sort -o $5
samtools view -F 4 -F 256 -F 2048 $5|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\t'$6'' > $7

# demo
# bowtie2-build nucl_nonerude.fna geneset_bowtie2
# bowtie2 -x geneset_bowtie2 -1 test.unmap.1.fastq.gz -2 test.unmap.2.fastq.gz -S bowtie2.sam -p 16
# samtools view bowtie2.sam -bS | samtools sort -o bowtie2.bam
## 过滤掉未比对到基因集的reads(Flag=4),过滤掉二次比对的reads(Flag=256),过滤掉嵌合reads(Flag=2048),过滤掉比对数小于2的
# samtools view -F 4 -F 256 -F 2048 bowtie2.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\tbowtie2' > bowtie2.count
