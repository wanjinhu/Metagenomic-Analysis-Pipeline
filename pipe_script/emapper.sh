#! /bin/bash

emapper.py -i $1 -o $2 --cpu 0 --usemem
cut -f1,12 $3|grep -v "^#"|sed 's/ko://g'|sed '1i gene\tko'|awk '$2 !~ /-/ {print}' > $4
cut -f1,13 $3|grep -v "^#"|sed '1i gene\tpathway'|awk '$2 !~ /-/ {print}' > $5

## demo
# 非冗余基因集prot_nonerude.faa利用eggnog-mapper功能注释, --cpu 0表示用所有cpu, --usemem表示断点续存
# emapper.py -i prot_nonerude.faa -o eggnog --cpu 0 --usemem
# 生成非冗余基因集的KO列表
# less -S eggnog.emapper.annotations|cut -f1,12|grep -v "^#"|sed 's/ko://g'|sed '1i gene\tko'|grep -v "-" > KEGG_KO.txt
# 生成非冗余基因集的PATHWAY列表
# less -S eggnog.emapper.annotations|cut -f1,13|grep -v "^#"|sed '1i gene\tpathway'|grep -v "-" > KEGG_PATHWAY.txt
# 生成非冗余基因集的COG列表
# less -S eggnog.emapper.annotations|grep -v '^#'|cut -f1,7|sed '1i gene\tCOG_code'|grep -v '-' > COG_CODE.txt
