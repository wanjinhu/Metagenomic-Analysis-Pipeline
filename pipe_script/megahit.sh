#! /bin/bash

megahit -1 $1 -2 $2 -o $3 --out-prefix $4 -t 8
seqkit seq -m 500 $5 --remove-gaps > $6
sed -i 's/>/>'$4'_/g' $6
