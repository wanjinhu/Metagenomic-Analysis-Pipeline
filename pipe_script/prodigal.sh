#! /bin/bash

prodigal -p meta -a $1 -m -d $2 -o $3 -f gff -s $4 -i $5
