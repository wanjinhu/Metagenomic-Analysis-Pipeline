#! /bin/bash

bwa mem -t 32 $1 $2 $3 | samtools view -bS - | samtools sort - > $4
# 过滤掉未比对到基因集的reads(Flag=4),过滤掉二次比对的reads(Flag=256),过滤掉嵌合reads(Flag=2048),过滤掉比对数小于2的
samtools view -F 4 -F 256 -F 2048 $4|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\t'$5'' > $6

# demo
## bwa构建索引
# bwa index nucl_nonerude.fna -p geneset_bwa
## 样本基因丰度计算
# bwa mem -t 4 total/geneset_bwa DBYX1FB.unmap.1.fastq.gz DBYX1FB.unmap.2.fastq.gz | samtools view -bS - | samtools sort - > DBYX1FB_mapping_geneset.bam
## 过滤掉未比对到基因集的reads(Flag=4),过滤掉二次比对的reads(Flag=256),过滤掉嵌合reads(Flag=2048),过滤掉比对数小于2的
# samtools view -F 4 -F 256 -F 2048 DBYX1FB_mapping_geneset.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\t"$2; else print $1"\t0"}'|sed '1i gene\tDBYX1FB' > DBYX1FB.count
## 不过滤比对数目小于2的
# samtools view -F 4 -F 256 -F 2048 DBYX1FB_mapping_geneset.bam|awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {print $1"\t"$2}'|sed '1i gene\tDBYX1FB' > DBYX1FB.count