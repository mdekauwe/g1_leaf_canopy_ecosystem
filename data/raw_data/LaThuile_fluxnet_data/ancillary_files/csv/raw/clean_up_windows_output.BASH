#!/bin/bash

for i in *.csv
do
    sed 's/^M//g' $i > tmp
    mv tmp $i
done

perl -p -e 's/^M//g' SummaryAnc_txt.csv > tmp ; mv tmp SummaryAnc_txt.csv 
