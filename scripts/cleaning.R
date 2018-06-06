#!/usr/bin/Rscript

# Cleaning DPRtab file to remove sequencing errors

# Input : .dprtab file
argv <- commandArgs(TRUE)
# argv <- c("regiella.DPRtab","regiella")
# setwd("~/vcf_parser")
out.name <- argv[2]

source("vcf_analysis_fun.r")
library(GenomicRanges,quietly = T)

### Loading input
dprtab <- read.table(argv[1],header=T,stringsAsFactors = F)
colnames(dprtab)[1:4] <- c("CHROM","POS","REF","ALT")

# Remove sites with 0 coverage
dprtab <- dprtab[apply(dprtab,1,function(line){return(sum(line=="0,0"))})==0,]


###
### Optionnal filtering of regions with homology with other symbionts
###
blastfiles <- list.files(paste0("~/ref_genome/blast/",out.name,"/"),full.names = T)

intervals.total <- GRanges()
for (file in blastfiles){
  cat(basename(file)," : ")
  intervals <- blast_to_grange(file)
  intervals.total <- append(intervals.total,intervals)
  dprtab <- filter_intervals(dprtab,intervals)
}
removed.length <- sum(width(ranges(reduce(intervals.total))))
cat("Genome length removed from analysis because of homology with other genomes : ",removed.length,"\n")

###
### Coverage filtering
###
cov.tab <- apply(dprtab[,5:ncol(dprtab)],2,function(col){return(rowSums(line2matrix(col)))})
cov.tab <- cov.tab[,colmin(cov.tab)!=colMeans(cov.tab)] # Remove outgroup added by vcf-merge
cov.tab <- cov.tab[rowmin(cov.tab)!=rowMeans(cov.tab),] # remove sites differing only from outgroup
median.cov <- apply(cov.tab,2,median)
cov.tab.norm <- sweep(cov.tab,2,median.cov,"/")

to.remove <- c(which(rowSums(cov.tab.norm<0.2)>0.75*ncol(cov.tab.norm)),  # Removing sites with >5 times median coverage or <5 times
               which(rowSums(cov.tab.norm>5)>0.75*ncol(cov.tab.norm)))

cat("Removed",length(to.remove),"sites because of very high/low coverage\n")
if (length(to.remove)>0){
  dprtab <- dprtab[-to.remove,]
}

### Filtering rare variants###############################

tab.list <- list()  # Allele counts for each site/sample

for (j in 5:ncol(dprtab)){
  tab <- line2matrix(dprtab[,j])
  ############# Filtering sequencing errors ############
  # Coverage filter
  tab[tab<=3] <- 0
  
  # Frequency filter
  freq <- rowmin(tab)/rowSums(tab)
  
  freq[is.na(freq)] <- 0
  tab[freq<0.1 & freq>0,] <- remove.min(tab[freq<0.1 & freq>0,])
  
  # Storing filtered coverages
  tab.list[[j-4]] <- tab
  ######################################################
}
names(tab.list) <- colnames(dprtab)[5:ncol(dprtab)]


# Removing samples with low coverage and high polymorphism

to.remove <- c()
for (j in 1:length(tab.list)){
  tab <- tab.list[[j]]
  nbsites <- sum(apply(tab,1,function(vec){return(vec[1]==vec[2])}))
  cat(paste0("sites with equal coverage for both alleles : ",nbsites,"\t"))
  if (nbsites/nrow(tab)>0.05){
    cat("Removing sample : ",names(tab.list)[j],"\n")
    to.remove <- c(to.remove,j)
  } else {cat("\n")}
}
if (length(to.remove)>0){
  tab.list <- tab.list[-to.remove]
}


### Output : cleaned dprtab

# Writing cleaned read counts
dprtab.cleaned <- cbind(dprtab[,1:4],do.call(cbind,lapply(tab.list,matrix2line)))
# Removing lines with ALT allele after cleaning
dprtab.cleaned <- dprtab.cleaned[Reduce("+",tab.list)[,2]!=0,]

write.table(dprtab.cleaned,file=paste0(argv[1],".cleaned"),quote = F,row.names = F)
