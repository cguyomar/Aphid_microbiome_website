#!/bin/bash

# $1 : fof location
# $2 : number of threads for each job

while read line
  do
    read file1 <<< $(echo $line | cut -d " " -f1)
    read matefile <<< $(echo $line | cut -d " " -f2)
    read ID <<< $(echo $line | cut -d " " -f3)
    qsub -pe make $2 -cwd ~/scripts/mapping_bwa_gz_multi.sh $file1 $matefile $ID $2
  done < $1

