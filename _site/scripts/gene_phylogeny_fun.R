
# Input : 
# contig name and fasta file of the reference genome
# Output : 
# Reference sequence for the requested contig

get_ref_seq <- function(contig,arg){
  fasta <- read.table(arg$ref)[,1]
  first <- which(fasta==paste0(">",contig))
  seq <- fasta[first:length(fasta)]
  last <- grep(pattern = ">",seq)[2]
  if (is.na(last)==F){ # Not last contig
    seq <- seq[1:(last-1)]
  }
  seq <- do.call(paste0,as.list(as.character(seq))[-1])
  seq <- unlist(strsplit(seq,split=""))
  return(seq)
}

# Input : 
# reference sequence, allele coverage, gene coordinates, and requested sample name
# Output : 
# Gene reference sequence altered to match sample genotype

get_alt_seq <- function(ref,dprtab,contig,start.gene,end.gene,sample_name){
  cat(paste0("Contig : ", contig,"\n Sample : ",sample_name))
  nb.change <- 0
  
  dprtab.contig <- dprtab[dprtab$CHROM==contig,]
  
  snps <- dprtab.contig[nchar(dprtab.contig$REF==1) & nchar(dprtab.contig$ALT)==1,]
  indels <- dprtab.contig[!(nchar(dprtab.contig$REF==1) & nchar(dprtab.contig$ALT)==1),]
  
  # Applying snps
  if (nrow(snps)>0){
    tab.snps <- line2matrix(snps[,colnames(snps)==sample_name])
    snps.alt <- snps[which(apply(tab.snps,1,which.max)==2),] # Snps to apply
    snps.alt <- snps.alt[snps.alt$POS>=start.gene & snps.alt$POS <=end.gene,]
    # cat(paste0("Number of SNPs added in gene: ",nrow(snps.alt),"\n"))
    nb.change <- nb.change+nrow(snps.alt)
    ref[snps.alt$POS] <- snps.alt$ALT
  }
  
  # Barcodes to track position of gene start and end
  ref[start.gene] <- paste(c(ref[start.gene],"START"),collapse = "_")
  ref[end.gene] <- paste(c(ref[end.gene],"END"),collapse = "_")
  original_length <- end.gene - start.gene
  tmp.start <- ref[start.gene]
  tmp.end <- ref[end.gene]
  
  # Dealing with indels
  if (dim(indels)[1]>0){
    tab.indels <- line2matrix(indels[,colnames(indels)==sample_name])
    indels.alt <- indels[which(apply(tab.indels,1,which.max)==2),] # Indels to apply
    indels.alt <- indels.alt[indels.alt$POS < end.gene+25 & indels.alt$POS > start.gene - 50,] # Only in gene neighborhood
    
    if (nrow(indels.alt)>0){
      nb.change <- nb.change+nrow(indels.alt)
      for (i in nrow(indels.alt):1){
        start <- indels.alt[i,]$POS
        end <- indels.alt[i,]$POS+nchar(indels.alt[i,]$REF)-1
        
        alt.seq <- strsplit(indels.alt[i,]$ALT,split="")[[1]]
        # Trim alt allele to insert
        alt.seq <- alt.seq[-c(1,length(alt.seq))]
        
        # We need to be careful if END or start are in indel
        test_end <- grep("END",ref[(start+1):(end-1)])
        if (length(test_end)>0 ){ #End of gene in indel
          if (is.na(alt.seq[test_end])){  # Deletion
            alt.seq[length(alt.seq)] <- paste0(alt.seq[length(alt.seq)],"_END")
          } else { # Insertion
            if (alt.seq[test_end]==substr(ref[start+test_end],1,1)){
              alt.seq[test_end] <- ref[start+test_end]
            } else{
              print("ERROR : end of gene in indel")
            }
          }
          
          test_start <- grep("START",ref[(start+1):(end-1)])
          if (length(test_start)>0 ){ #End of gene in indel
            if (alt.seq[length(alt.seq) - (end-start-1) +test_start]==substr(ref[start+test_start],1,1)){
              alt.seq[length(alt.seq) - (end-start-1) +test_start] <- ref[start+test_start]
            } else{
              print("ERROR : start of gene in indel")
            }
          }
          ref <- c(ref[1:(start)],
                   alt.seq,
                   ref[(end):length(ref)])
          
          if (length(which(ref==tmp.start))==0 | length(which(ref==tmp.end))==0 ){
            print("### Warning : Start or end of gene in indel ###")
          }
          start.gene <- grep("_START",ref)
          end.gene <- grep("_END",ref)
        }
        # cat(paste0("Total gene length change : ", end.gene-start.gene- original_length,"\n \n"))
      }
    }
  }
  cat(paste0("   ",nb.change," variants \n"))
  ref[grepl(pattern = "START",ref)] <- strsplit(ref[grepl(pattern = "START",ref)],"_")[[1]][1]
  ref[grepl(pattern = "END",ref)] <- strsplit(ref[grepl(pattern = "END",ref)],"_")[[1]][1]
  
  return(paste0(ref[start.gene:(end.gene)],collapse = ""))
}


# Each Ipg file can contain several lines (with identical genomic coordinates)
# This func returns only one, preferably RefSeq sequences
filter_tab <- function(gene,strain){
  # Return only one ID if several provided
  tab <- read.csv(gene,sep="\t",header=T)
  tab <- tab[tab$Strain==strain,]
  if (nrow(tab)>1){
    if ("RefSeq" %in% tab$Source){
      tab <- tab[tab$Source=="RefSeq",]
      if (nrow(tab)>1) {
        if (length(unique(tab$Start))==1 & length(unique(tab$Stop))==1 ){
          tab <- tab[1,]
        } else {
          print(tab)
        }
      }
    }
  }
  return(tab)
}

rowmin <- function(tab){return(apply(tab,1,min))}

line2matrix <- function(line){
  return(matrix(as.numeric(unlist(strsplit(as.character(line),","))),ncol=2,byrow = T))
}

