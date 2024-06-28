rm(list = ls(all=T))

#set working directory
setwd("~/Desktop/Work/mLOY_Exomes/")


library(tidyverse)
library("Biostrings")
library("MatchIt")



## read the file saved from previous step (1_choose_chrY_exons.R) with selected exons on chrY
chYexons <- read.table("resources/xgen_plus_spikein_chrY_singleCopyGenes.GRCh38.bed", sep="\t")
colnames(chYexons) <- c("CHR", "start", "end")
## read exome capture (here the one used in the UKB)

## To selects control exons we will use only autosomal genes that rarely become parts of mosaic Chromosomal abberations
## As identified in Loh et al. Nature, 2018. The data with mCAs are available as individual level data
## as part of the UKB, andd available to approved UKB researchers.

## We used these data to calculate for each exon in the exome capture how often it becomes
## part of different types of mCAs - gains, loss, loss of heterozygosity, or unknown.
## The result is saved to resources/xgen_plus_spikein.GRCh38.mcas.count.bed

#open this file

exome <- read.table("resources/xgen_plus_spikein.GRCh38.mcas.count.bed", sep="\t", header=T)

## We leave the exons that have the count of gains and losses (copy change) less or equal to 10 (arbitrary)

exome <- exome %>% filter(!(CHR %in% c("chrX", "chrY"))) %>% filter((n_gain <= 10) & (n_loss <= 10))



##########
## Now for each selected chrY exon we will select 100 autosomal exons matched on  GC content and length

## Calculate length

exome$length <- exome$end - exome$start
chYexons$length <- chYexons$end - chYexons$start


## Calculate GC content
## For that you need human hg38 genome fasta file 

genome <- readDNAStringSet("resources/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta")

# for ausosomal
GC_cont <- c()
for (i in 1:dim(exome)[1]) {
  chr=exome$CHR[i]
  start <- exome$start[i]
  end <- exome$end[i]
  
  
  piece <- subseq(genome[grepl(paste0(chr,"  AC"), names(genome))], start = start, end = end)
  gc_content <- letterFrequency(piece, letters="GC")/letterFrequency(piece, letters="ATGC")
  GC_cont <- c(GC_cont, gc_content)}

exome$GC <- GC_cont
rm(piece, GC_cont)

# for chrY

GC_cont <- c()
for (i in 1:dim(chYexons)[1]) {
  chr=chYexons$CHR[i]
  start <- chYexons$start[i]
  end <- chYexons$end[i]
  
  
  piece <- subseq(genome[grepl(paste0(chr,"  AC"), names(genome))], start = start, end = end)
  gc_content <- letterFrequency(piece, letters="GC")/letterFrequency(piece, letters="ATGC")
  GC_cont <- c(GC_cont, gc_content)}

chYexons$GC <- GC_cont
rm(piece, GC_cont, genome)


# merge 
exome <- rbind(exome[,c(1,2,3,8,9)], chYexons)


## Match
exome$chrY <- ifelse(exome$CHR == "chrY", 1, 0)
exome$GC <- round(exome$GC, digits = 2)

matched <- matchit(chrY ~ GC+length, data = exome,
        method = "nearest", caliper = c(GC=0.01, length=5), ratio=100, discard="control")

## Look at the matching summary

#summary(matched)

### Save the result
## This file already exists for the ukb


write.table(match.data(matched)[,c(1:3,9)], "resources/exome_capture_matched_intervals.groups.bed", row.names = F, col.names = F, sep="\t", quote = F)
write.table(match.data(matched)[,c(1:3)], "resources/exome_capture_matched_intervals.bed", row.names = F, col.names = F, sep="\t", quote = F)







