#!/bin/bash

source /local/env/envR.sh

echo Symbiote > left
cut -f1 -d " " ArPo28.symb.comptage >> left
#echo unmapped >> left
while read file
	do echo $file > $file.symb.temp
	cut -d " " -f2 $file.symb.comptage >> $file.symb.temp
done < samples.list
paste -d ";" left $(sed ':a;N;$!ba;s/\n/.symb.temp /g' samples.list ) > symb.comptages
rm left
rm *.temp

Rscript sort_tab.R