#!/bin/bash

source /local/env/envsamtools-1.6.sh
source /local/env/envjava.sh

filename=$(basename $1)
outname=$(echo $filename | cut -d "." -f1)

samtools view -b -f 12 -F 256 $1 > ./${outname}_unmapped.bam

samtools bam2fq -1 ${outname}_unmapped.1.fastq -2 ${outname}_unmapped.2.fastq ${outname}_unmapped.bam

java -jar ~/soft/Trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 1 -phred33 ${outname}_unmapped.1.fastq ${outname}_unmapped.2.fastq ${outname}_unmapped_cleaned.1.fastq ${outname}_unmapped_unpaired.1.fastq ${outname}_unmapped_cleaned.2.fastq ${outname}_unmapped_unpaired.2.fastq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
