#! /bin/bash

grep -E '(p__)|(clade_name)' $1 |grep -v 'c__'|sed 's/|/;/g' > $2
grep -E '(c__)|(clade_name)' $1 |grep -v 'o__'|sed 's/|/;/g' > $3
grep -E '(o__)|(clade_name)' $1 |grep -v 'f__'|sed 's/|/;/g' > $4
grep -E '(f__)|(clade_name)' $1 |grep -v 'g__'|sed 's/|/;/g' > $5
grep -E '(g__)|(clade_name)' $1 |grep -v 's__'|sed 's/|/;/g' > $6
grep -E '(s__)|(clade_name)' $1 |grep -v 't__'|sed 's/|/;/g' > $7
