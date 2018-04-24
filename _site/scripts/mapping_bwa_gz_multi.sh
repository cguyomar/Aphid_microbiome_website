#!/bin/bash

source /local/env/envbwa-0.7.10.sh
source /local/env/envsamtools.sh

# $1 and $2 : paired end libraries
# $3 : outname
# $4 number of threads

echo "processed file : " $3
bwa mem -t $4 /omaha-beach/cguyomar/ref_genome/ref_genome_review $1 $2 | samtools view -b - > $3.bam

samtools sort -@ $4 $3.bam $3_sorted
mv $3_sorted.bam $3.bam

samtools index $3.bam
samtools flagstat $3.bam > $3.flagstat

