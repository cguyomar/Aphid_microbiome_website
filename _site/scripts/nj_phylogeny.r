#!/usr/bin/Rscript

# Analysis of variant calling results

# Input : .dprtab file
argv <- commandArgs(TRUE)
# argv <- "regiella"
# setwd("~/vcf_parser/")
name <- argv[1]

source("vcf_analysis_fun.r")
library(ggplot2)
library(reshape2)
library(gplots)
library(RColorBrewer)
library(dendextend)
library(GenomicRanges)
library(phytools)

# Colours for each biotype
cols <- c("#6593be","#7e3f34","#cb4472","#c6a349","#cb4ebf","#d55936","#7bac65","#d290a5","#4c5d2f","#8b7ed4","#70c64b","#5f3a70","#673fbe","#5db7a9")

###
### Loading input
###

dprtab <- read.table(paste0(name,".DPRtab.cleaned"),header=T,stringsAsFactors = F)

###
tab.list <- list()  # Allele counts for each site/sample
maj.tab <- c()      # Major allele for each site/sample
freq.tab <- c()     # minor allele frequency

for (j in 5:ncol(dprtab)){
  tab <- line2matrix(dprtab[,j])
  
  freq <- rowmin(tab)/rowSums(tab)
  freq[is.na(freq)] <- 0

  tab.list[[j-4]] <- tab
}
names(tab.list) <- colnames(dprtab)[5:ncol(dprtab)]



###
### Computing genotype for each sample from major alleles ###############
###


for (j in 1:length(tab.list)){
  tab <- tab.list[[j]]

  freq <- rowmin(tab)/rowSums(tab)
  freq[is.na(freq)] <- 0
  freq.tab <- cbind(freq.tab,freq)

  maj <- ifelse(apply(tab,1,which.max)==1,"REF","ALT")
  # maj[tab[,1]==tab[,2]] <- sample(c("REF","ALT"),size = sum(tab[,1]==tab[,2]),replace = T) # Random assignation
  maj.tab <- cbind(maj.tab,maj)
}

colnames(freq.tab) <- colnames(maj.tab) <- names(tab.list)
rownames(freq.tab) <- rownames(maj.tab) <- paste(1:nrow(freq.tab),dprtab[,1],dprtab[,2],sep="_")




###
### Computing distance matrix for nj phylogeny ##################
###

allele.mat <- maj.tab

# Change names
tab <- read.table("~/vcf_parser/change_names.csv",sep=",",row.names=1)
# allele.mat <- allele.mat[,colnames(allele.mat) %in% rownames(tab)]
colnames(allele.mat)[colnames(allele.mat) %in% rownames(tab)] <- as.character(tab[colnames(allele.mat)[colnames(allele.mat) %in% rownames(tab)] ,])



### Nj inference
dist.mat <- get_distmat(allele.mat)
nj.tree <- nj(dist.mat)
nj.tree <- midpoint.root(nj.tree)

# output statistics 
cat("Number of variants (including outgroup) :",nrow(maj.tab),"\n")
cat("Max distance (including outgroup) :",max(dist.mat),"\n")

tmp <- maj.tab[,-ncol(maj.tab)]
tmp <- maj.tab[,-c(47,52)]
tmp <- tmp[apply(tmp,1,function(line){return(length(unique(line)))})>1,]
cat("Number of variants (without outgroup) : ",nrow(tmp),"\n")
cat("Max distance (excluding outgroup) :",max(dist.mat[-nrow(dist.mat),-nrow(dist.mat)]),"\n")




plot(nj.tree)

boot.values <- boot.phylo(nj.tree,x=t(allele.mat),FUN = function(xx) {midpoint.root(nj(get_distmat(t(xx))))})
nj.tree$node.label <- boot.values

svg(file=paste0(name,"_njtree.svg"),width=12,height = 12)
plotBreakLongEdges(nj.tree,tip.color = tips_colors(nj.tree$tip.label),align.tip.label = T,
                   cex=1.3,edge.width=3,n=2)
nodelabels(nj.tree$node.label,bg=ifelse(as.numeric(nj.tree$node.label)>90,"yellowgreen","brown2"))
add.scale.bar()

write.tree(nj.tree,file=paste0("results/",name,"_nj.newick"))
save(dist.mat,nj.tree,boot.values,file=paste0("results/",name,".Rdata"))
