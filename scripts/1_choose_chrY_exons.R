### Choose chrY exons to calculate the normalized chrY coverage

#load libraries
library(tidyverse)

#set working directory
setwd("~/Desktop/Work/mLOY_Exomes/")

#exome caprute intervals (.bed)
#example for exome capture used in the UKB
#replace with the relevant intervals file
exome_capture <- read.table("resources/ukb_xgen_plus_spikein.GRCh38.bed")
exome_capture <- exome_capture[exome_capture$V1 == "chrY",]

#Homo sapiens genome annotation
annot <- read.table(gzfile("resources/Homo_sapiens.GRCh38.112.chr.gtf.gz"), header=F, sep='\t')
annot <- annot[(annot$V1 == "Y") & (annot$V3 == "gene"),]
annot <- annot %>%
  mutate(V9 = strsplit(annot$V9, "; ")) %>%
  unnest(V9) %>%
  group_by(V4) %>%
  mutate(row = row_number()) %>%
  spread(row, V9)

colnames(annot) <- c("Chr", 'source', 'type', 'start', 'end', 'V6', 'strand', 'V8', 'EnsemblID',
                     'gene_version', 'gene_name', 'gene_source', 'gene_biotype')

annot <- annot[!grepl("gene_biotype" ,annot$gene_source),]


annot <- annot %>%
  mutate(EnsemblID = gsub("gene_id ", "", EnsemblID)) %>%
  mutate(gene_version = gsub("gene_version ", "", gene_version)) %>%
  mutate(gene_name =gsub("gene_name ", "", gene_name)) %>%
  mutate(gene_source = gsub("gene_source ", "", gene_source)) %>%
  mutate(gene_biotype = gsub("gene_biotype |;", "", gene_biotype))

# we are interested only in protein coding genes
annot <- annot[annot$gene_biotype == "protein_coding",]
#and only in the ones located in MSY
annot <- annot[(annot$start > 2781479) & (annot$end > 2781479) & (annot$start < 56887903) &
                 (annot$end < 56887903),]


# The list of single-copy genes (based on Skaletsky et al., 2003, https://doi.org/10.1038/nature01722)

single_genes <- c("SRY" ,    "RPS4Y1",  "ZFY" ,"AMELY" ,  "TBL1Y" , "PRKY",  "USP9Y" ,
                  "DDX3Y" ,  "UTY"  ,   "TMSB4Y" , "NLGN4Y",  "KDM5D" ,  "EIF1AY" , "RPS4Y2" )

annot_filtered <- annot[annot$gene_name %in% single_genes,]

#Filter exons only in those 14 single-copy genes

single_genes_coord <- c()

for (i in 1:nrow(annot_filtered)) {
  single_genes_coord <- c(single_genes_coord, annot_filtered[i,]$start:annot_filtered[i,]$end)
}

exome_capture_single_genes <- exome_capture[exome_capture$V2 %in% single_genes_coord & exome_capture$V3 %in% single_genes_coord,]


#chosen <- read.table("resources/xgen_plus_spikein_chrY_singleCopyGenes_annotated.GRCh38.bed")
#exome_capture_single_genes <- merge(exome_capture_single_genes, chosen, by=c("V1", "V2", "V3"), all.x=T)
#write.table(exome_capture_single_genes, "resources/test.annot.bed")

#write the result
write.table(exome_capture_single_genes, "resources/xgen_plus_spikein_chrY_singleCopyGenes.GRCh38.bed", quote=F, sep="\t")

## In this example with exome capture used in the ukb we end up with 184 exons. 
## In our paper we further filtered this list of exons by looking at the coverage of each exon
## in people who have no evidence of mLOY as identified by PAR-LOY method. We excluded exons that on average
## deviated higher than median +- igr.
## The final list of  chosenexons on chrY is saved in resources/exome_capture_matched_intervals.bed
## (together with matched autosomal exons). It can be directry used starting from step3 (3_run_mosdepth.sh) in case the exome capture 
##is the same as in the UKB. But the method should perform well enough without this additional optimization step.







