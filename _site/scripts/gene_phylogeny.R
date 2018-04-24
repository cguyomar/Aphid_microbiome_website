#!/usr/bin/Rscript

##  Retrive genomic corrdinates from a protein refseq ID


### Libraries and functions 
library(IRanges)
library(bio3d)
source("./gene_phylogeny_fun.R")

argv <- commandArgs(TRUE)
# argv <- c("regiella")
# setwd("~/vcf_parser/")

name <- argv[1]



# Sets of parameters -------------------------------------------------

args.list <- list()
args.list$buchnera <- list(dir="~/protein_fetching/buchnera/",
                           strain="APS",
                           ref="~/ref_genome/genomes/buchnera.fasta",
                           contig="Buchnera_gi15616630refNC_002528.1")
args.list$hamiltonella <- list(dir="~/protein_fetching/hamiltonella/",
                               strain="T5A",
                               ref="~/ref_genome/genomes/hamiltonella.fasta",
                               contig="Hamiltonella_gi238897251refNC_012751.1")
args.list$regiella <- list(dir="~/protein_fetching/regiella/",
                           strain="5.15",
                           ref="~/ref_genome/genomes/regiella.fasta")
args.list$serratia <- list(dir="~/protein_fetching/serratia/",
                           strain="Tucson",
                           ref="~/ref_genome/genomes/serratia.fasta")


###
### Loading input
###
print(paste0(name,".DPRtab.cleaned"))
dprtab <- read.table(paste0(name,".DPRtab.cleaned"),header=T,stringsAsFactors = F)
head(dprtab)

eval(parse(text=paste0("arg <- args.list$",name)))

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



###################   Removing overlapping Indels #########################

# We need to parse dprtab to remove overlapping indels (and keep only one)
for (contig in unique(dprtab$CHROM)){
  # print(contig)
  # Subsetting for a contig
  which.indels <- which(dprtab$CHROM==contig & (nchar(dprtab$REF)!=1 | nchar(dprtab$ALT)!=1))
  indels <- dprtab[which.indels,]
  indels$END <- indels$POS+nchar(indels$REF)-1
  
  # Finding overlaps
  indels_range <- IRanges(indels$POS,indels$END)
  ol <- findOverlaps(indels_range,ignoreSelf=T,ignoreRedundant=F)
  
  # Resolving overlaps
  if (length(ol)>0){
    to.remove <- c()
    while (length(ol)>0){
      to.remove <- c(to.remove,which.indels[which.max(as.table(ol))])
      which.indels <- which.indels[-which.max(as.table(ol))]
      indels <- indels[-which.max(as.table(ol)),]
      indels_range <- IRanges(indels$POS,indels$END)
      ol <- findOverlaps(indels_range,ignoreSelf=T,ignoreRedundant=F)
    }
    
    # Applying to dprtab
    dprtab <- dprtab[-to.remove,]
  }
}




directory <- arg$dir

fasta <- read.table(arg$ref)
strain <- arg$strain

genes <- list.files(directory,pattern = "*.res")
genes.tab <- matrix()
setwd(directory)

tab <- do.call(rbind,lapply(as.list(genes),filter_tab,strain))



###
### Computing haplotypes
###

aln.list <- list()

for (i in 1:nrow(tab)){
  print(i)
  gene <- tab[i,]
  start.gene <- as.numeric(gene$Start)
  end.gene <- as.numeric(gene$Stop)
  
  if (is.null(arg$contig)){ # Several contigs genome
    contig <- as.character(gene$Nucleotide.Accession)
    line.contig <- which(sapply(fasta,function(line) {return(grepl(gsub("_","",contig),gsub("_","",line)))}))
    contig <- substring(as.character(fasta[line.contig,]),2,nchar(as.character(fasta[line.contig,])))
  } else {
    contig <- arg$contig
  }
  
  ref <- get_ref_seq(contig,arg)
  for (j in 5:ncol(dprtab)){
    sample_name <- colnames(dprtab)[j]
    print(sample_name)

    seq.alt <- get_alt_seq(ref,dprtab,contig,start.gene,end.gene,sample_name)
    seq.alt <- strsplit(seq.alt,split="")[[1]]
    
    if (j == 5){ # First sample
      aln <- seqbind(ref[start.gene:end.gene],seq.alt)
    } else {
      aln <- seqbind(aln,seq.alt)
    }
  }
  aln.list <- c(aln.list,list(aln))
}

###
### Write each gene in a fasta file
###
for (i in 1:length(aln.list)){
  fasta <- as.fasta(aln.list[[i]]$ali,id = c("REF",colnames(dprtab)[5:ncol(dprtab)]))
  bio3d::write.fasta(fasta,file=paste0("alignements/",tab$Protein[i],".fasta"))
}




