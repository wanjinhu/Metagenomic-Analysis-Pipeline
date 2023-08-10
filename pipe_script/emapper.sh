#! /bin/bash

emapper.py -i $1 -o $2 --cpu 0 --usemem
cut -f1,12 $3|grep -v "^#"|sed 's/ko://g'|sed '1i gene\tko'|awk '$2 !~ /-/ {print}' > $4
cut -f1,13 $3|grep -v "^#"|sed '1i gene\tpathway'|awk '$2 !~ /-/ {print}' > $5
