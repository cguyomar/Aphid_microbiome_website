#!/bin/bash

source /local/env/envsamtools.sh

filename=$(basename $1)
outname=$(echo $filename | cut -d "." -f1)

samtools view -b $1 APSE1_ENA\|AF157835\|AF157835.1 > /omaha-beach/cguyomar/genomecov/APSE1/$outname.APSE1.bam
samtools view -b $1 Buchnera_gi15616630refNC_002528.1 > /omaha-beach/cguyomar/genomecov/buchnera/$outname.buchnera.bam
samtools view -b $1 Hamiltonella_gi238897251refNC_012751.1> /omaha-beach/cguyomar/genomecov/hamiltonella/$outname.hamiltonella.bam
samtools view -b $1 Mitochondrion_gi213948225refNC_011594.1 > /omaha-beach/cguyomar/genomecov/mitochondrie/$outname1.mito.bam
samtools view -b $1 $(cat /omaha-beach/cguyomar/genomecov/regiella/regiella.contigs) > /omaha-beach/cguyomar/genomecov/regiella/$outname.regiella.bam 
samtools view -b $1 $(cat /omaha-beach/cguyomar/ref_genome/genomes/rickettsia.contigs) > /omaha-beach/cguyomar/genomecov/rickettsia/$outname.rickettsia.
bam
samtools view -b $1 Rickettsiella_viridis > /omaha-beach/cguyomar/genomecov/rickettsiella/$outname.rickettsiella.bam 
samtools view -b $1 $(cat /omaha-beach/cguyomar/genomecov/serratia/serratia.contigs) > /omaha-beach/cguyomar/genomecov/serratia/$outname.serratia.bam 
samtools view -b $1 $(cat /omaha-beach/cguyomar/genomecov/spiroplasma/spiro.contigs) > /omaha-beach/cguyomar/genomecov/spiroplasma/$outname.spiroplasma.
bam 
samtools view -b $1 $(cat /omaha-beach/cguyomar/ref_genome/genomes/fukatsuia.contigs) > /omaha-beach/cguyomar/genomecov/fukatsuia/$outname.fukatsuia.bam
samtools view -b $1.bam Mitochondrion_gi213948225refNC_011594.1 > /omaha-beach/cguyomar/genomecov/mitochondrie/$1.mito.bam
