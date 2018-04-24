#!/bin/bash

source /local/env/envsamtools.sh

filename=$(basename $1)
outname=$(echo $filename | cut -d "." -f1)

samtools idxstats $1 | sort -k 1b,1  > $outname.idxstats
join $outname.idxstats /omaha-beach/cguyomar/scaffolds_review.list > $outname.scaff.comptage
awk '{SUM+=$4} END {print "unmapped "SUM}' $outname.idxstats >> $outname.scaff.comptage
awk '{a[$5]+=$3}END{for (i in a) print i,a[i]}' $outname.scaff.comptage | sed 1d > $outname.symb.comptage
awk '{SUM+=$4} END {print "unmapped "SUM}' $outname.idxstats >> $outname.symb.comptage

