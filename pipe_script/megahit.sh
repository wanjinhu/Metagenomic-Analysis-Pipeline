#! /bin/bash

megahit -1 $1 -2 $2 -o $3 --out-prefix $4 -t 32
seqkit seq -m 500 $5 --remove-gaps > $6
sed -i 's/>/>'$4'_/g' $6
## demo
# megahit组装
# megahit -1 DBYX1FB.unmap.1.fastq.gz -2 DBYX1FB.unmap.2.fastq.gz -o DBYX1FB_megahit --out-prefix DBYX1FB -t 8
# 过滤500bp以下的组装序列,修改序列名称
# seqkit seq -m 500 DBYX1FB_megahit/DBYX1FB.contigs.fa --remove-gaps > DBYX1FB.contigs_500.fa
# sed -i 's/>/>DBYX1FB_/g' DBYX1FB.contigs_500.fa