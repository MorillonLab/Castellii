#  ***********
#  *  TOOLS  *
#  ***********

#  R version 3.4.4
#  bedtools v2.27.0

#  ***************
#  *  Functions  * 
#  ***************

source("functions.R")

#  **************
#  *  libraries *
#  **************

library(reshape2)
library(ggplot2)
library(VennDiagram)
library(DESeq2)

#  *****************
#  *  ANNOTATIONS  *
#  *****************

### N. castellii annotation from https://fungi.ensembl.org

annot <- readGff("Naumovozyma_castellii_cbs_4309.ASM23734v1.38.Chr.gff3")
annot <- annot[annot$Type %in% c("gene", "snRNA", "snoRNA", "tRNA", "rRNA"),]

# annotation of aslncRNA from Alcid & Tsukiyama. Nat Struct Mol Biol. 2016 May; 23(5): 450–455. 

aslncRNA <- readGff("castellii_aslncRNA.gff")

# annotation of new lncRNA detected in this study

DUT_file <- "DE_new_transcripts_xrn1_dcr1/transcrits_319/new_DUT.gff"
DUT <- readGff(DUT_file)

SUT_file <- "DE_new_transcripts_xrn1_dcr1/transcrits_319/new_SUT.gff"
SUT <- readGff(SUT_file)

XUT_file <- "DE_new_transcripts_xrn1_dcr1/transcrits_319/new_XUT.gff"
XUT <- readGff(XUT_file)

CUT_file <- "DE_new_transcripts_rrp6/transcrits_159/new_CUT.gff"  
CUT <- readGff(CUT_file)

# Refining SUTs & CUTs (Materials and Methods : Annotation of lncRNAs) :
# 116 CUT overlapped >50% of transcripts defined as SUTs, these 116 transcripts were not considered as SUTs, but CUTs.

# intersect SUT and CUT (using bedtools)

system(paste0("bedtools intersect -nonamecheck -s -wo -a ", CUT_file,
              "-b ", SUT_file," > intersect_new_CUT_vs_SUT.tab"))

intersect_CUT_SUT <- readIntersect("intersect_new_CUT_vs_SUT.tab")

# 50% overlap between SUT and CUT

intersect_CUT_SUT_50 <- intersect_CUT_SUT[intersect_CUT_SUT$V19/(intersect_CUT_SUT$V5 - intersect_CUT_SUT$V4 + 1) >= 0.5 
                                          & intersect_CUT_SUT$V19/(intersect_CUT_SUT$V14 - intersect_CUT_SUT$V13 + 1) >= 0.5,]

# SUT removed

SUT <- SUT[!(SUT$ID %in% as.character(intersect_CUT_SUT_50$V18)),]

writeGff(df = SUT, out = "SUT_corrected.gff")

SUT_file <- "SUT_corrected.gff"

# merge all annotations

annot_all <- rbind.data.frame(annot, XUT, CUT, SUT, aslncRNA, DUT)
annot_all <- annot_all[order(annot_all$Chr, annot_all$Start, annot_all$Stop),]
writeGff(df = annot_all, out = "annot_all.gff")



### S. cerevisiae annotations

# annotation from http://www.yeastgenome.org/

annot_cer <- readGff("saccharomyces_cerevisiae.gff")
annot_cer <- annot_cer[annot_cer$Type %in% c("gene", "ncRNA", "snoRNA", "snRNA", "tRNA", "rRNA"),]

# XUT annotation van Dijk et al. Nature. 2011 Jun 22;475(7354):114-7.

XUT_cer <- readGff("XUT_20150527_corrected_by_CAGEseq.gff")

# CUT and SUT annotation from Xu et al. Nature. 2009 Feb 19; 457(7232): 1033–1037. 

tmp <- readGff("NUTs_blasted_from_Steinmetz.gff")
CUT_cer <- tmp[tmp$Type == "CUT",]
SUT_cer <- tmp[tmp$Type == "SUT",]

# merge annotations
annot_all_cer <- rbind.data.frame(annot_cer, XUT_cer, CUT_cer, SUT_cer)
annot_all_cer <- annot_all_cer[order(annot_all_cer$Chr, annot_all_cer$Start, annot_all_cer$Stop),]
writeGff(df = annot_all_cer, out = "annot_all_cer.gff")


#  *************************
#  *  ANALYSIS ANNOTATION  *
#  *************************

### N. castellii antisens non-coding transcripts

# SUT antisens

system(paste0("bedtools intersect -nonamecheck -S -wo -a ", SUT_file,
                "-b annot_all.gff > SUT_antisens_all_annot.tab"))
SUT_antisens_all_annot <- readIntersect("SUT_antisens_all_annot.tab")
SUT_antisens_mRNA <- SUT_antisens_all_annot[SUT_antisens_all_annot$V12 == "gene",]

# CUT antisens

system(paste0("bedtools intersect -nonamecheck -S -wo -a ", CUT_file,
              "-b annot_all.gff > CUT_antisens_all_annot.tab"))

CUT_antisens_all_annot <- readIntersect("CUT_antisens_all_annot.tab")
CUT_antisens_mRNA <- CUT_antisens_all_annot[CUT_antisens_all_annot$V12 == "gene",]

# XUT antisens

system(paste0("bedtools intersect -nonamecheck -S -wo -a ", XUT_file,
              "-b annot_all.gff > XUT_antisens_all_annot.tab"))

XUT_antisens_all_annot <- readIntersect("XUT_antisens_all_annot.tab")
XUT_antisens_mRNA <- XUT_antisens_all_annot[XUT_antisens_all_annot$V12 == "gene",]

# DUT antisens

system(paste0("bedtools intersect -nonamecheck -S -wo -a ", DUT_file,
              "-b annot_all.gff > DUT_antisens_all_annot.tab"))

DUT_antisens_all_annot <- readIntersect("DUT_antisens_all_annot.tab")
DUT_antisens_mRNA <- DUT_antisens_all_annot[DUT_antisens_all_annot$V12 == "gene",]

# all aslncRNA antisens of coding transcripts

nc.anti.1nt <- c(unique(SUT_antisens_mRNA[SUT_antisens_mRNA$V19 >= 1,]$V9), 
                 unique(CUT_antisens_mRNA[CUT_antisens_mRNA$V19 >= 1,]$V9), 
                 unique(XUT_antisens_mRNA[XUT_antisens_mRNA$V19 >= 1,]$V9),
                 unique(DUT_antisens_mRNA[DUT_antisens_mRNA$V19 >= 1,]$V9),)

### overlap XUT and CUT 

system(paste0("bedtools intersect -nonamecheck -s -wo -a ", XUT_file, "-b ", CUT_file," > overlap_XUT_CUT.tab"))

overlap_XUT_CUT <- readIntersect("overlap_XUT_CUT.tab")
overlap_XUT_CUT_50 <- overlap_XUT_CUT[overlap_XUT_CUT$V19/(overlap_XUT_CUT$V5 - overlap_XUT_CUT$V4 + 1) >= 0.5 & 
                                        overlap_XUT_CUT$V19/(overlap_XUT_CUT$V14 - overlap_XUT_CUT$V13 + 1) >= 0.5,]

# overlap new non-coding RNA detected and previously (Alcid & Tsukiyama)

system(paste0("cat ", XUT_file, " ", SUT_file," ", CUT_file," | ",
              "bedtools sort -i stdin | ",
              "bedtools intersect -nonamecheck -s -wo -a stdin -b castellii_aslncRNA.gff > ",
              "overlap_nc_aslnc_alcid.tab"))
overlap_nc_aslnc_alcid <- readIntersect("overlap_nc_aslnc_alcid.tab")


### S. cerevisiae antisens non-coding transcripts

# SUT cerevisiae antisens

system(paste0("cat annot_all_cer.gff | awk 'OFS=\"\t\"{if($3 == \"SUT\") print $0}' | ",
              "bedtools intersect -nonamecheck -S -wo -a stdin ",
              "-b annot_all_cer.gff > SUT_antisens_all_annot_cer.tab"))

SUT_antisens_all_annot_cer <- readIntersect(file = "SUT_antisens_all_annot_cer.tab")
SUT_antisens_mRNA_cer <- SUT_antisens_all_annot_cer[SUT_antisens_all_annot_cer$V12 == "gene",]

# CUT cerevisiae antisens
system(paste0("cat annot_all_cer.gff | awk 'OFS=\"\t\"{if($3 == \"CUT\") print $0}' | ",
              "bedtools intersect -nonamecheck -S -wo -a stdin ",
              "-b annot_all_cer.gff > CUT_antisens_all_annot_cer.tab"))

CUT_antisens_all_annot_cer <- readIntersect(file = "CUT_antisens_all_annot_cer.tab")
CUT_antisens_mRNA_cer <- CUT_antisens_all_annot_cer[CUT_antisens_all_annot_cer$V12 == "gene",]

# XUT cerevisiae antisens
system(paste0("cat annot_all_cer.gff | awk 'OFS=\"\t\"{if($3 == \"XUT\") print $0}' | ",
              "bedtools intersect -nonamecheck -S -wo -a stdin ",
              "-b annot_all_cer.gff > XUT_antisens_all_annot_cer.tab"))

XUT_antisens_all_annot_cer <- readIntersect(file = "XUT_antisens_all_annot_cer.tab")
XUT_antisens_mRNA_cer <- XUT_antisens_all_annot_cer[XUT_antisens_all_annot_cer$V12 == "gene",]


#  ******************************
#  *  COUNT TABLE FROM RNA-SEQ  *
#  ******************************

### make count table for WT, xrn1, dcr1, xrn1 dcr1 total RNA-seq data (this study)

# bam files

file_names = ... # vector with path to bam files here

# experimental design : sample names

description_data = c("WT_1", "WT_2", "xrn1_1", "xrn1_2", 
                     "dcr1_1", "dcr1_2", "xrn1_dcr1_1", "xrn1_dcr1_2")

# raw bamHandler

datarough=readBam(file_names,ncore=8,libraryType="inverse")

# ERCC counts for normalization
# number of reads mapping on ERCC spike-in used

gff_ERCC = readGff("ERCC.gff")
gff_ERCC$ID = gff_ERCC$Chr

counts_ERCC = countReadsBams(data = datarough, description_data = description_data, annot = gff_ERCC, ncore = ncore)

poids_xrn1 = colSums(counts_ERCC[,grepl(".*_[1-2]\ readcount", colnames(counts_ERCC))]) 

# normalized bamHandler
# normalization coefficient formula : min/poids 

data=normalizeBams(datarough, poids_xrn1)

# for each condtion, take samples average
# WT
data[[9]]=meanBam(data[[1]], data[[2]])
# xrn1
data[[10]]=meanBam(data[[3]], data[[4]])
# dcr1
data[[11]]=meanBam(data[[5]], data[[6]]) 
# xrn1 dcr1
data[[12]]=meanBam(data[[7]], data[[8]])

# pooled sample names
description_data[9:12]= c("WT", "xrn1", "dcr1", "xrn1_dcr1") 

# raw counts
counts_xrn1_dcr1_raw <- countReadsBams(data = data, description_data = description_data,
                                       annot = annot_all, ncore=4, normalized=FALSE)
colnames(counts_xrn1_dcr1_raw) <- make.names(colnames(counts_xrn1_dcr1_raw)) 
counts_xrn1_dcr1_raw <- cbind.data.frame(annot_all, counts_xrn1_dcr1_raw)

# normalized counts
counts_xrn1_dcr1 <- countReadsBams(data = data, description_data = description_data,
                                   annot = annot_all, ncore=4, normalized=TRUE)
colnames(counts_xrn1_dcr1) <- make.names(colnames(counts_xrn1_dcr1)) 
counts_xrn1_dcr1 <- cbind.data.frame(annot_all, counts_xrn1_dcr1)


### make count table for WT, rrp6, dcr1, rrp6 dcr1 total RNA-seq data (Alcid & Tsukiyama)
# bam files

file_names = ... # vector with path to bam files here

# experimental design : sample names

description_data = c("WT_1", "WT_2", "rrp6_1", "rrp6_2", 
                     "dcr_1", "dcr_2", "dcr1_rrp6_1", "dcr1_rrp6_2")

# raw bamHandler

datarough=readBam(file_names,ncore=ncore,libraryType="inverse")

# counts on coding transcripts for normalization
# number of reads mapping on coding transcripts used

counts_coding <- countReadsBams(data = datarough, 
                                description_data = description_data[1:8], 
                                annot = annot[annot$Type == "gene",], 
                                ncore = ncore, 
                                normalized = F)
poids_rrp6 <- colSums(counts_coding[,grep("readcount", colnames(counts_coding))])

# normalized bamHandler

data=normalizeBams(datarough, poids_rrp6)

# replicate mean
# WT
data[[9]]=meanBam(data[[1]], data[[2]])
# xrn1
data[[10]]=meanBam(data[[3]], data[[4]])
# dcr1
data[[11]]=meanBam(data[[5]], data[[6]]) 
# xrn1 dcr1
data[[12]]=meanBam(data[[7]], data[[8]])

# pooled sample names

description_data[9:12]= c("WT", "rrp6", "dcr1", "dcr1_rrp6") 

# raw counts

counts_rrp6_dcr1_raw <- countReadsBams(data = data, description_data = description_data,
                                       annot = annot_all, ncore=4, normalized=FALSE)
colnames(counts_rrp6_dcr1_raw) <- make.names(colnames(counts_rrp6_dcr1_raw)) 
counts_rrp6_dcr1_raw <- cbind.data.frame(annot_all, counts_rrp6_dcr1)

# normalized counts

counts_rrp6_dcr1 <- countReadsBams(data = data, description_data = description_data,
                                   annot = annot_all, ncore=4, normalized=TRUE)
colnames(counts_rrp6_dcr1) <- make.names(colnames(counts_rrp6_dcr1)) 
counts_rrp6_dcr1 <- cbind.data.frame(annot_all, counts_rrp6_dcr1)


### make count table for small RNA-seq (this study)

# bam file (22-23nt mapping on N. castellii)

files.Ncas <- ... # vector with path to bam files here

# sample names

description_data <- c("WT_1", "WT_2", "xrn1_1", "xrn1_2", "dcr1_1", "dcr1_2", "xrn1_dcr1_1", "xrn1_dcr1_2", "Dcr1_GFP")

# make bamHandler

datarough <- readBam(file = files.Ncas, libraryType = "inverse", ncore = 8)

# counts raw

counts_small_raw <- countReadsBams(data = datarough, 
                                   description_data = description_data, 
                                   annot = annot_all, 
                                   ncore = 8, normalized = F)
colnames(counts_small_raw) <- make.names(colnames(counts_small_raw))
counts_small_raw <- cbind.data.frame(annot_all, 
                                     counts_small_raw)

# mapping stats, obtained from script processing_small_RNAseq.sh

mapping_stats <- read.csv("mapping_stats_smallRNA.tab"), 
                          header = T, sep = "\t", stringsAsFactors = F)

# coeff pombe

nReads.Spom.cen <- mapping_stats$nReads.Spom.cen.22.23
coeff.pom <- min(nReads.Spom.cen)/nReads.Spom.cen

# normalized counts
counts_small <- matrix(unlist(lapply(description_data, 
                                     function(sample){
                                       i_coeff <- coeff.pom[which(description_data == sample)]
                                       readcount <- counts_small_raw[,paste0(sample, ".readcount")]*i_coeff
                                       densities <- readcount/(annot_all$Stop - annot_all$Start + 1)
                                       return(list(readcount, densities))
                                     })), nrow = nrow(annot_all), ncol = length(description_data)*2)
colnames(counts_small) <- paste0(rep(description_data, each = 2), c(".readcount", ".densities"))
counts_small <- cbind.data.frame(annot_all, counts_small)


#  ******************
#  *  PLOT FIGURES  *
#  ******************

### Figure 1.B, 1.C, 1.D & 1.E : Heatmap and density plot of the expression fold-change of (coding and non-coding) transcripts relative to the corresponding WT strain (total RNA-seq)

# compute fold-change

df <- cbind.data.frame(ID = as.character(counts_xrn1_dcr1$ID),
                       Type = counts_xrn1_dcr1$Type,
                       log2_fc_dcr_wt = log2((counts_xrn1_dcr1$dcr1.readcount+1)/
                                              (counts_xrn1_dcr1$WT.readcount+1)),
                       log2_fc_xrn1_wt = log2((counts_xrn1_dcr1$xrn1.readcount+1)/
                                               (counts_xrn1_dcr1$WT.readcount+1)),
                       log2_fc_xrn1_dcr1_wt = log2((counts_xrn1_dcr1$xrn1_dcr1.readcount+1)/
                                                      (counts_xrn1_dcr1$WT.readcount+1)),
                       log2_fc_rrp6_wt = log2((counts_rrp6_dcr1$rrp6.readcount+1)/
                                                (counts_rrp6_dcr1$WT.readcount+1)),
                       log2_fc_rrp6_dcr1_wt = log2((counts_rrp6_dcr1$dcr1_rrp6.readcount+1)/
                                                      (counts_rrp6_dcr1$WT.readcount+1)),
                       stringsAsFactors = F)
  
# is lncRNA antisens of ORF ?

is_anti <- ifelse(df$ID %in% nc.anti.1nt,
                  T, F)
df <- cbind.data.frame(df, is_anti)
  
# Fig 1.B : Heatmap lncRNA

mat <- c()
types <- c("XUT", "CUT", "SUT", "DUT")
  
# make fold-change matrix, where lncRNA are grouped by type (CUT, SUT, XUT, DUT, solo, antisens)

for(t in types){
  if(t %in% c("XUT", "CUT", "SUT")){
    for(is_anti_ in c(T, F)){
      mat_ <- df[df$Type == t & df$is_anti == is_anti_, ]
      rownames(mat_) <- mat_$ID
      mat_ <- mat_[,c(1,2,3,4,6,8)]
      mat_ <- mat_[rowSums(mat_[,3:6]) != 0,]
      clust <- hclust(as.dist(1 - cor(t(mat_[,3:6]), use = "pa", method = "pearson")), 
                      method = "ward.D2")
      mat_ <- mat_[clust$order,]
        
      if(is.null(mat)){
        mat <- mat_
      } else {
        mat <- rbind(mat,
                     mat_)
      }
    }
  } else {
    mat_ <- rbind(df[df$Type == t & df$is_anti == F,],
                  df[df$Type == t & df$is_anti == T,])
    rownames(mat_) <- mat_$ID
    mat_ <- mat_[,c(1,2,3,4,6,8)]
    mat_ <- mat_[rowSums(mat_[,3:6]) != 0,]
    
    if(is.null(mat)){
      mat <- mat_
    } else {
      mat <- rbind(mat,
                   mat_)
    }
  }
}
  
# for colorscale

x = quantile(seq(-10,7,1), probs = seq(0,1, 1/length(-10:7)))
ncol = length(x)
  
# plot heatmap

layout(mat = matrix(data = 1:7, nrow = 7, ncol = 1), 
       heights = c(0.05, 0.07, 0.16, 0.27, 0.15, 0.24, 0.06))
par(mar = c(0.5,4,0.5,4))

# SUT

image(x = 1:(ncol(mat[mat$Type == "SUT" & mat$is_anti == F,3:5])), 
      y = 1:nrow(mat[mat$Type == "SUT" & mat$is_anti == F,3:5]), 
      z = as.matrix(t(mat[mat$Type == "SUT" & mat$is_anti == F,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "s", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "SUT" & mat$is_anti == F,]), side = 4, line = 1, las = 2)
  
image(x = 1:(ncol(mat[mat$Type == "SUT" & mat$is_anti == T,3:5])), 
      y = 1:nrow(mat[mat$Type == "SUT" & mat$is_anti == T,3:5]), 
      z = as.matrix(t(mat[mat$Type == "SUT" & mat$is_anti == T,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "as", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "SUT" & mat$is_anti == T,]), side = 4, line = 1, las = 2)

# CUT

image(x = 1:(ncol(mat[mat$Type == "CUT" & mat$is_anti == F,3:5])), 
      y = 1:nrow(mat[mat$Type == "CUT" & mat$is_anti == F,3:5]), 
      z = as.matrix(t(mat[mat$Type == "CUT" & mat$is_anti == F,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "s", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "CUT" & mat$is_anti == F,]), side = 4, line = 1, las = 2)
  
image(x = 1:(ncol(mat[mat$Type == "CUT" & mat$is_anti == T,3:5])), 
      y = 1:nrow(mat[mat$Type == "CUT" & mat$is_anti == T,3:5]), 
      z = as.matrix(t(mat[mat$Type == "CUT" & mat$is_anti == T,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "as", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "CUT" & mat$is_anti == T,]), side = 4, line = 1, las = 2)

# XUT

image(x = 1:(ncol(mat[mat$Type == "XUT" & mat$is_anti == F,3:5])), 
      y = 1:nrow(mat[mat$Type == "XUT" & mat$is_anti == F,3:5]), 
      z = as.matrix(t(mat[mat$Type == "XUT" & mat$is_anti == F,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "s", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "XUT" & mat$is_anti == F,]), side = 4, line = 1, las = 2)

image(x = 1:(ncol(mat[mat$Type == "XUT" & mat$is_anti == T,3:5])), 
      y = 1:nrow(mat[mat$Type == "XUT" & mat$is_anti == T,3:5]), 
      z = as.matrix(t(mat[mat$Type == "XUT" & mat$is_anti == T,3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
box()
abline(v = 1.5)
abline(v = 2.5)
mtext(text = "as", side = 2, las = 2, line = 1)
mtext(text = nrow(mat[mat$Type == "XUT" & mat$is_anti == T,]), side = 4, line = 1, las = 2)

# DUT

image(x = 1:(ncol(mat[mat$Type == "DUT",3:5])), 
      y = 1:nrow(mat[mat$Type == "DUT",3:5]), 
      z = as.matrix(t(mat[mat$Type == "DUT",3:5])), 
      col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      breaks = x, 
      xlab = "", ylab = "",
      xaxt = "n", yaxt = "n")
abline(v = 1.5)
abline(v = 2.5)
abline(h = 4.5)
box()
mtext(text = "s", side = 2, line = 1, at = 9.5, las = 2)
mtext(text = "as", side = 2, line = 1, at = 4.5, las = 2)
mtext(text = nrow(mat[mat$Type == "DUT" & mat$is_anti == F,]), side = 4, line = 1, las = 2, at = 9.5)
mtext(text = nrow(mat[mat$Type == "DUT" & mat$is_anti == T,]), side = 4, line = 1, las = 2, at = 4.5)

# heatmap legend

image(x = x, y = 1:2,
      z = as.matrix(x), col = colorRampPalette(c("navyblue","white", "red4"), bias = 0.8)(ncol-1),
      xlab = "log2 fold-change mutant/wt", ylab = "", yaxt = "n")
  

### Fig 1.C, 1.D, 1.E : Density plot

dens_list <- list()
  
# for color and layout

cols <- c("cornflowerblue", "black", "chocolate4", 
          rgb(255,0,0, maxColorValue = 255), 
          rgb(0,176,80, maxColorValue = 255), 
          rgb(127,127,127, maxColorValue = 255))
names(cols) <- c("gene", "snoRNA", "tRNA", "XUT", "CUT", "SUT")
  
xlims <- rbind(c(-3, 7),
                 c(-4, 7),
                 c(-2,2.5))
rownames(xlims) <- c("xrn1", "rrp6", "dcr")
  
# density plot for each mutant (xrn1, rrp6, dcr1) vs WT 

for(fc in c("xrn1", "rrp6", "dcr")){
  dens_list[[fc]] <- list()
  tmp_x <- c()
  tmp_y <- c()
  for(t in c("gene", "snoRNA", "tRNA", "XUT", "CUT", "SUT")){
    dens_list[[fc]][[t]] <- density(df[df$Type == t, paste0("log2_fc_",fc,"_wt")])
    tmp_x <- c(tmp_x, dens_list[[fc]][[t]]$x)
    tmp_y <- c(tmp_y, dens_list[[fc]][[t]]$y)
  }
  
  x_lims <- xlims[fc,]
  y_lims <- c(0, max(tmp_y))
    
  plot(x = x_lims, y = y_lims, col = "white",
       xlab = paste0(fc,"/WT ratio (log2)"), ylab = "Density", bty = "n")
  for(t in c("gene", "snoRNA", "tRNA", "XUT", "CUT", "SUT")){
    lines(x = dens_list[[fc]][[t]]$x, y = dens_list[[fc]][[t]]$y,
          col = cols[t], lwd = 1.5)
  }
  abline(h = 0)
  legend("topleft", legend = c("ORF", "snoRNA", "tRNA", "XUT", "CUT", "SUT"),
         cols, bty = "n", text.col = cols)
}


### Fig 1.F, S1.H & 2.A : 
# Box-plot of densities for antisense and solo non-coding transcripts in WT cells (total RNA-seq)
# Box-plot of densities for antisense and solo XUT in xrn1 cells (total RNA-seq)
# Box-plot of densities for the antisense and solo XUTs in the xrn1 and xrn1 dcr1 strains (total RNA-seq)

# make data frame with lncRNA densities for all strains

df <- counts_xrn1_dcr1[counts_xrn1_dcr1$Type %in% c("SUT","XUT", "CUT"),
                       c(1:6,grep("WT.densities|xrn1.densities|dcr1.densities|xrn1.dcr1.densities", 
                                  colnames(counts_xrn1_dcr1)))]

df <- cbind.data.frame(df, is.anti = df$ID %in% nc.anti.1nt)
  
df <- melt(data = df, measure.vars = 7:10, variable.name = "strain", value.name = "densities")
  
df$strain <- gsub(".densities", "", df$strain)
df$strain <- factor(df$strain, levels = c("WT", "dcr1", "xrn1", "xrn1_dcr1"))
  
df$Type <- factor(df$Type, levels = c("SUT", "CUT", "XUT"))
  
df$is.anti <- factor(df$is.anti, levels = c(TRUE, FALSE))
 
# Fig 1.F : solo and antisens SUT, CUT & XUT in WT cells

p <- ggplot(data = df[df$strain %in% c("WT") & 
                        df$Type %in% c("SUT","XUT","CUT"),], 
            mapping = aes(x = Type, y = log2(densities), fill = is.anti)) + 
  geom_boxplot(lwd = 1, outlier.colour = NA, notch = T) +
  scale_fill_manual(values = c("grey80", "grey40"), labels = c("antisens", "solo"), name = "") + 
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.3), 
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "RNA seq expression of non-coding transcripts in WT") + 
  labs(x = "", y = "RNA-seq, log2 densities")
p 

# Fig S1.H : antisense and solo XUTs densities in xrn1 cells (total RNA-seq)

p <- ggplot(data = df[df$strain %in% c("xrn1") & 
                        df$Type == "XUT",], 
            mapping = aes(x = is.anti, y = log2(densities), fill = strain)) + 
  geom_boxplot(lwd = 1, outlier.colour = NA, notch = T) +
  scale_fill_manual(values = c("red")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.3), 
         plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("antisens", "solo"))+
  ggtitle(label = "RNA seq expression of XUT") + 
  labs(x = "", y = "RNA-seq, log2 densities")
p

# Fig 2.A : antisense and solo XUTs in the xrn1 and xrn1 dcr1 strains

p <- ggplot(data = df[df$strain %in% c("xrn1", "xrn1_dcr1"),], 
            mapping = aes(x = is.anti, y = log2(densities), fill = strain)) + 
  geom_boxplot(lwd = 1, outlier.colour = NA, notch = T) +
  scale_fill_manual(values = c("grey80", "grey40")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.3), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("antisens", "solo"))+
  ggtitle(label = "RNA seq expression of XUT") + 
  labs(x = "", y = "RNA-seq, log2 densities")
p  

# Pairwise Wilcoxon test for each non-coding type (solo and antisens) in all strains

for(t in c("SUT", "CUT", "XUT")){
  for.test <- df[df$Type == t,]
  for.test$densities <- log2(for.test$densities)
  for.test <- for.test[is.finite(for.test$densities),]
  w_test <- pairwise.wilcox.test(x = for.test$densities, 
                                 g = paste0(for.test$strain, "_", for.test$is.anti),
                                 p.adjust.method = "BH")
  table <- format(x = w_test$p.value, digits = 3)
  colnames(table) <- gsub("_FALSE", " solo", colnames(table))
  colnames(table) <- gsub("_TRUE", " antisens", colnames(table))
  rownames(table) <- gsub("_FALSE", " solo", rownames(table))
  rownames(table) <- gsub("_TRUE", " antisens", rownames(table))
    
  f <- paste0("Wilcoxon_pairwise_comparison_",t,"_xrn1_dcr1.csv")
  
  cat("#", w_test$method, "for ", t,"\n", file = f)
  cat("#Adjust method : Benjamini & Hochberg\n", file = f, append = T)
  cat("\n", file = f, append = T)
  cat("\t", file = f, append = T)
  write.table(table, file = f, 
              col.names = T, row.names = T, na = "-", quote = F, sep = "," , append = T)
}

### Fig S1.A, S1.B, S1.C & S1.D : Scatter plot of tag density coding and non-coding transcripts (total RNA-seq)
# WT vs dcr1, WT vs xrn1, WT vs rrp6 & dcr1 vs xrn1 

# make list with transcripts densities

count_list <- list()
count_list[[1]] <- counts_xrn1_dcr1[,c("Type", "WT.densities", "xrn1.densities")]
count_list[[2]] <- counts_rrp6_dcr1[,c("Type", "WT.densities", "rrp6.densities")]
count_list[[3]] <- counts_xrn1_dcr1[,c("Type", "WT.densities", "dcr1.densities")]
count_list[[4]] <- counts_xrn1_dcr1[,c("Type", "dcr1.densities", "xrn1.densities")]
  
cols <- c("grey80", "black",  
          "red", 
          rgb(0,176,80, maxColorValue = 255), 
          rgb(127,127,127, maxColorValue = 255),
          "red")
names(cols) <- c("gene", "snoRNA", "XUT", "CUT", "SUT","DUT")

bgs <- cols
names(bgs) <- c("gene", "snoRNA", "XUT", "CUT", "SUT","DUT")
  
type_list <- list()
type_list[[1]] <- c("gene", "snoRNA", "SUT", "XUT")
type_list[[2]] <- c("gene", "snoRNA", "SUT", "CUT")
type_list[[3]] <- c("gene", "snoRNA", "SUT", "DUT")
type_list[[4]] <- c("gene", "DUT")
  
for(i in 1:length(count_list)){
    
  s1 <- gsub(".densities", "", colnames(count_list[[i]])[2])
  s2 <- gsub(".densities", "", colnames(count_list[[i]])[3])
    
  plot(x = log2(count_list[[i]][count_list[[i]]$Type %in% names(cols),2]),
       y = log2(count_list[[i]][count_list[[i]]$Type %in% names(cols),3]), col = "white", 
       xlab = paste0(s1, ", density (tag/nt, log2)"),
       ylab = paste0(s2, ", density (tag/nt, log2)"), bty = "n")
  for(t in type_list[[i]]){
    points(x = log2(count_list[[i]][count_list[[i]]$Type == t,2]),
           y = log2(count_list[[i]][count_list[[i]]$Type == t, 3]),
           bg = bgs[t], col = cols[t], pch = 21, cex = 0.7)
  }

  abline(0,1, col = "red")
  legend("topleft", legend = c("ORF", type_list[[i]][-1]), text.col = cols[type_list[[i]]], bty = "n")  
}


### Fig S1.E : Venn diagram showing the overlap (≥ 50%) between CUTs and XUTs.

grid.newpage()
draw.pairwise.venn(area1 = nrow(XUT), area2 = nrow(CUT), cross.area = nrow(overlap_XUT_CUT_50), 
                   category = c("XUT", "CUT"))


### Fig S1.F : Overlap (≥ 1 nt) between the SUTs, XUTs and CUTs identified in this work and the 170 previously annotated aslncRNAs (Alcid and Tsukiyama 2016)

mat <- cbind(c(length(unique(overlap_nc_aslnc_alcid[overlap_nc_aslnc_alcid$V3 == "SUT","V9"])), nrow(SUT)),
             c(length(unique(overlap_nc_aslnc_alcid[overlap_nc_aslnc_alcid$V3 == "XUT", "V9"])), nrow(XUT)),
             c(length(unique(overlap_nc_aslnc_alcid[overlap_nc_aslnc_alcid$V3 == "CUT", "V9"])), nrow(CUT)))
barplot(mat, beside = T, col = c("grey80", "grey40"), 
        names.arg = c("SUTs", "XUTs", "CUTs"), ylab = "# transcripts",
        main = "overlap with previously annotated aslncRNA")
legend("topleft", legend = c("overlap previously annotated", "new"), bty = "n", 
       fill = c("grey80", "grey40"))


### Fig S1.G : Proportion of SUTs, CUTs and XUTs that are antisense (≥ 1 nt) to an ORF or to any annotated transcript 

mat <- cbind(c(length(unique(SUT_antisens_all_annot[SUT_antisens_all_annot$V12 == "gene","V9"]))/nrow(SUT), 
               length(unique(SUT_antisens_all_annot[,"V9"]))/nrow(SUT)),
             c(length(unique(CUT_antisens_all_annot[CUT_antisens_all_annot$V12 == "gene","V9"]))/nrow(CUT), 
               length(unique(CUT_antisens_all_annot[,"V9"]))/nrow(CUT)),
             c(length(unique(XUT_antisens_all_annot[XUT_antisens_all_annot$V12 == "gene","V9"]))/nrow(XUT), 
               length(unique(XUT_antisens_all_annot[,"V9"]))/nrow(XUT)))

barplot(mat, beside = T, col = c("white", "black"), ylim = c(0,1), 
        ylab = "Proportion", names.arg = c("SUT", "CUT", "XUT"), 
        main = "Proportion of aslncRNA")


### Fig 2.B & S1.I : 
# Box-plot of densities for the antisense and solo CUTs in the rrp6 and rrp6 dcr1 strains (total RNA-seq)
# Box-plot of densities for the antisense and solo CUTs in rrp6 cells (total RNA-seq)

df <- counts_rrp6_dcr1[counts_rrp6_dcr1$Type %in% c("CUT"),
                        c(1:6,grep("WT.densities|rrp6.densities|dcr1.densities|dcr1_rrp6.densities", 
                                   colnames(counts_rrp6_dcr1)))]

df <- cbind.data.frame(df, is.anti = df$ID %in% nc.anti.1nt)

df <- melt(data = df, measure.vars = 7:10, variable.name = "strain", value.name = "densities")

df$strain <- gsub(".densities", "", df$strain)
df$strain <- factor(df$strain, levels = c("WT", "dcr1", "rrp6", "dcr1_rrp6"))
  
df$Type <- factor(df$Type, levels = c("CUT"))
  
df$is.anti <- factor(df$is.anti, levels = c(TRUE, FALSE))
  
# Fig 2.B : antisense and solo CUTs in the rrp6 and rrp6 dcr1 strains

p <- ggplot(data = df[df$strain %in% c("rrp6", "dcr1_rrp6"),], 
            mapping = aes(x = is.anti, y = log2(densities), fill = strain)) + 
  geom_boxplot(lwd = 1, outlier.colour = NA, notch = T) +
  scale_fill_manual(values = c("grey80", "grey40")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.3), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("antisens", "solo"))+
  ggtitle(label = "RNA seq expression of CUT") + 
  labs(x = "", y = "RNA-seq, log2 densities") + 
  lims(y = c(-6, 0))
p  
  
# Fig S1.I : antisense and solo CUTs in rrp6 cells

p <- ggplot(data = df[df$strain %in% c("rrp6") & 
                        df$Type == "CUT",], 
            mapping = aes(x = is.anti, y = log2(densities), fill = strain)) + 
  geom_boxplot(lwd = 1, outlier.colour = NA, notch = T) +
  scale_fill_manual(values = c("green")) + 
  theme(panel.background = element_blank(),
        axis.line = element_line(size = 0.3), 
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("antisens", "solo"))+
  ggtitle(label = "RNA seq expression of CUT") + 
  labs(x = "", y = "RNA-seq, log2 densities") + 
  lims(y = c(-5, 0))
p

# Wilcoxon test

for(t in c("CUT")){
  for.test$densities <- log2(for.test$densities)
  for.test <- for.test[is.finite(for.test$densities),]
  for.test$cond <- paste0(for.test$strain, "_", for.test$is.anti)
  for.test$cond <- factor(for.test$cond, levels = unique(for.test$cond))
  w_test <- pairwise.wilcox.test(x = for.test$densities, 
                                 g = for.test$cond,
                                 p.adjust.method = "BH")
  table <- format(x = w_test$p.value, digits = 3)
  colnames(table) <- gsub("_FALSE", " solo", colnames(table))
  colnames(table) <- gsub("_TRUE", " antisens", colnames(table))
  rownames(table) <- gsub("_FALSE", " solo", rownames(table))
  rownames(table) <- gsub("_TRUE", " antisens", rownames(table))
   
  f <- paste0("Wilcoxon_pairwise_comparison_",t,"_rrp6_dcr1.csv")
  
  cat("#", w_test$method, "for ", t,"\n", file = f)
  cat("#Adjust method : Benjamini & Hochberg\n", file = f, append = T)
  cat("\n", file = f, append = T)
  cat("\t", file = f, append = T)
  write.table(table, file = f, 
              col.names = T, row.names = T, na = "-", quote = F, sep = "," , append = T)
}


### Fig S2.E : Number of Dcr1-sensitive protein-coding genes in the presence or absence of Xrn1 or Rrp6 (total RNA-seq)

# Differential expression analysis for rrp6 vs rrp6 dcr strains
  
outDE <- paste0("DE_analysis_rrp6_dcr1_vs_rrp6/")
dir.create(outDE)

# prepare data for DESeq2

sample_data <- cbind.data.frame(Names = grep("_[12].readcount", 
                                             colnames(counts_rrp6_dcr1_raw), value = T),
                                Condition = gsub("_[12].readcount", "", grep("_[12].readcount", 
                                                                             colnames(counts_rrp6_dcr1_raw), value = T)))
    
counts_for_dds <- counts_rrp6_dcr1_raw[counts_rrp6_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT") ,
                                        as.character(sample_data$Names)]
rownames(counts_for_dds) <- counts_rrp6_dcr1_raw[counts_rrp6_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT"),]$ID

dds <- DESeqDataSetFromMatrix(countData = counts_for_dds,
                              colData = sample_data,
                              design = ~Condition)
    
# use normalisation factors computed

sizeFactors(dds) <- 1/(min(poids_rrp6)/poids_rrp6)

# run DESeq2

dds <- DESeq(dds)
    
# get results

res <- results(object = dds, contrast = c("Condition", "dcr1_rrp6", "rrp6"), alpha = 0.05)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
    
counts_ <- counts_rrp6_dcr1[counts_rrp6_dcr1$ID %in% rownames(res),]

# compute fold-change

fc <- (counts_$dcr1_rrp6.readcount+1) / (counts_$rrp6.readcount+1) 
    
# get table up, down and unchanged transcripts

table_up <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,
                                     c("ID", grep("rrp6.*densities", colnames(counts_), value = T))], 
                             fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2],
                             padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,]$padj)
write.table(table_up, file = paste0(outDE, "table_", t,"_up_rrp6_dcr1_vs_rrp6.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
table_down <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,
                                       c("ID", grep("rrp6.*densities", colnames(counts_), value = T))], 
                               fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5],
                               padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,]$padj)
write.table(table_down, file = paste0(outDE, "table_gene_down_rrp6_dcr1_vs_rrp6.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
table_unch <- cbind.data.frame(counts_[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),
                                       c("ID", grep("rrp6.*densities", colnames(counts_), value = T))], 
                               fc = fc[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5))],
                               padj = res[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),]$padj)
write.table(table_unch, file = paste0(outDE, "table_gene_unchanged_rrp6_dcr1_vs_rrp6.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)

# get number of genes affected
n_genes_de_rrp6_dcr1_vs_rrp6 <- nrow(table_down) + nrow(table_up)


# Differential expression analysis for xrn1 vs xrn1 dcr strains
 
outDE <- paste0("DE_analysis_xrn1_dcr1_vs_xrn1/")
dir.create(outDE)
 
sample_data <- cbind.data.frame(Names = grep("_[12].readcount", 
                                             colnames(counts_xrn1_dcr1_raw), value = T),
                                Condition = gsub("_[12].readcount", "", grep("_[12].readcount", 
                                                                             colnames(counts_xrn1_dcr1_raw), value = T)))
  
counts_for_dds <- counts_xrn1_dcr1_raw[counts_xrn1_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT") ,
                                        as.character(sample_data$Names)]
rownames(counts_for_dds) <- counts_xrn1_dcr1_raw[counts_xrn1_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT"),]$ID
    
dds <- DESeqDataSetFromMatrix(countData = counts_for_dds,
                                  colData = sample_data,
                                  design = ~Condition)
    
sizeFactors(dds) <- 1/(min(poids_xrn1)/poids_xrn1)

dds <- DESeq(dds)
    
res <- results(object = dds, contrast = c("Condition", "xrn1_dcr1", "xrn1"), alpha = 0.05)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
    
counts_ <- counts_xrn1_dcr1[counts_xrn1_dcr1$ID %in% rownames(res),]
    
fc <- (counts_$xrn1_dcr1.readcount+1) / (counts_$xrn1.readcount+1) 
    
table_up <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,
                                     c("ID", grep("xrn1.*densities", colnames(counts_), value = T))],
                             fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2],
                             padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,]$padj)
write.table(table_up, file = paste0(outDE, "table_gene_up_xrn1_dcr1_vs_xrn1.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
table_down <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,
                                       c("ID", grep("xrn1.*densities", colnames(counts_), value = T))], 
                               fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5],
                               padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,]$padj)
write.table(table_down, file = paste0(outDE, "table_gene_down_xrn1_dcr1_vs_xrn1.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
table_unch <- cbind.data.frame(counts_[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),
                                       c("ID", grep("xrn1.*densities", colnames(counts_), value = T))], 
                               fc = fc[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5))],
                               padj = res[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),]$padj)
write.table(table_unch, file = paste0(outDE, "table_gene_unchanged_xrn1_dcr1_vs_xrn1.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
n_genes_de_xrn1_dcr1_vs_xrn1 <- nrow(table_down) + nrow(table_up)

  
# Differential expression analysis for WT vs dcr strains
 
outDE <- paste0("DE_analysis_dcr1_vs_WT/")
dir.create(outDE)
    
sample_data <- cbind.data.frame(Names = grep("_[12].readcount", 
                                             colnames(counts_xrn1_dcr1_raw), value = T),
                                Condition = gsub("_[12].readcount", "", grep("_[12].readcount",
                                                                             colnames(counts_xrn1_dcr1_raw), value = T)))
    
counts_for_dds <- counts_xrn1_dcr1_raw[counts_xrn1_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT") ,
                                        as.character(sample_data$Names)]
rownames(counts_for_dds) <- counts_xrn1_dcr1_raw[counts_xrn1_dcr1_raw$Type %in% c("gene", "CUT", "SUT", "XUT","DUT"),]$ID

dds <- DESeqDataSetFromMatrix(countData = counts_for_dds,
                              colData = sample_data,
                              design = ~Condition)
    
sizeFactors(dds) <- 1/(min(poids_xrn1)/poids_xrn1)
dds <- DESeq(dds)
    
res <- results(object = dds, contrast = c("Condition", "dcr1", "WT"), alpha = 0.05)
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
    
counts_ <- counts_xrn1_dcr1[counts_xrn1_dcr1$ID %in% rownames(res),]
    
fc <- (counts_$dcr1.readcount+1) / (counts_$WT.readcount+1) 
    
table_up <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,
                                     c("ID", grep("WT.densities|^dcr1.densities", colnames(counts_), value = T))], 
                             fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2],
                             padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc >= 2,]$padj)
write.table(table_up, file = paste0(outDE, "table_gene_up_dcr1_vs_WT.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
   
table_down <- cbind.data.frame(counts_[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,
                                       c("ID", grep("WT.densities|^dcr1.densities", colnames(counts_), value = T))], 
                               fc = fc[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5],
                               padj = res[res$padj < 0.05 & counts_$Type == "gene" & fc <= 0.5,]$padj)
write.table(table_down, file = paste0(outDE, "table_gene_down_dcr1_vs_WT.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
table_unch <- cbind.data.frame(counts_[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),
                                       c("ID", grep("WT.densities|^dcr1.densities", colnames(counts_), value = T))], 
                               fc = fc[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5))],
                               padj = res[counts_$Type == "gene" & (res$padj >= 0.05 | (fc < 2 & fc > 0.5)),]$padj)
write.table(table_unch, file = paste0(outDE, "table_gene_unchanged_dcr1_vs_WT.tab"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
    
n_genes_de_dcr1_vs_WT <- nrow(table_down) + nrow(table_up) 
  
# barplot number of Dcr1-sensitive protein-coding genes in the presence or absence of Xrn1 or Rrp6

barplot(c(n_genes_de_dcr1_vs_WT,  n_genes_de_xrn1_dcr1_vs_xrn1,  n_genes_de_rrp6_dcr1_vs_rrp6), 
        names.arg = c("dcr vs WT", "xrn1 dcr vs xrn1", "rrp6 dcr vs rrp6"), col = "grey60")


### Fig.3 A : WT/dcr1 FC for solo and antisens SUT, CUT and XUT (small RNA)

# data

df <- cbind.data.frame(annot_all, 
                       log2_fc = log2((rowMeans(counts_small[,grep("WT.*readcount", colnames(counts_small))]) + 1) / 
                                        (rowMeans(counts_small[,grep("^dcr.*readcount", colnames(counts_small))]) + 1)))
df <- df[df$Type %in% c("SUT", "CUT", "XUT"),]
df$is_anti <- df$ID %in% nc.anti.1nt
df$Type <- factor(df$Type, levels = c("SUT", "CUT", "XUT"))

# boxplot

par(mar = c(5,4,4,2))
cols <- c(rgb(127,127,127, maxColorValue = 255), rgb(0,176,80, maxColorValue = 255), rgb(255,0,0, maxColorValue = 255))
boxplot(log2_fc ~ Type + is_anti, data = df, 
        outline = F, frame = F, xaxt = "n", notch = T, 
        at = c(1,2,3,5,6,7), col = cols, 
        lwd = 2, medlwd = 3, staplewex = 0, lty = 1, 
        ylab = "small RNA-seq, log2 fold-change (WT/dcr1)")
axis(side = 1, at = c(2,6), labels = c("solo", "antisens"), tick = F)
legend("topleft", legend = c("SUTs", "CUTs", "XUTs"), 
       text.col = cols, text.width = 3, bty = "n")

# Wilcoxon test

for(anti in c(T, F)){
  
  w_test <- pairwise.wilcox.test(x = df[df$is_anti == anti,]$log2_fc, 
                                 g = df[df$is_anti == anti,]$Type)
  table <- format(x = w_test$p.value, digits = 3)
  
  
  f <- paste0("Wilcoxon_pairwise_comparison_small_RNA_nc_",
              ifelse(anti, "antisens", "solo"),
              "_nc_fc_wt_over_dcr1.csv")
  
  cat("#", w_test$method,"\n", file = f)
  cat("# Small RNA expression fold-change (WT/dcr1, log2) vs ",
      ifelse(anti, "antisens", "solo"),
      " non-coding type (SUT, CUT, XUT)\n", sep = " ", file = f, append = T)
  cat("# Adjust method : Benjamini & Hochberg\n", file = f, append = T)
  cat("\n", file = f, append = T)
  cat("\t", file = f, append = T)
  write.table(table, file = f, 
              col.names = T, row.names = T, na = "-", quote = F, sep = "," , append = T)
}


### Fig S3.A : size distribution and 1st base (small RNAseq)

# size distribution and 1st base
# file *_first_base_small_no_tRNA_rRNA.tab were obtained from script processing_small_RNAseq.sh

# make list of reformated tables with size of read & first base

list_table_size_distribution <- list()

for(index_ in index){
  
  t <- read.csv(file = paste0( index_, "_first_base_small_no_tRNA_rRNA.tab"), 
                header = F, sep = "\t", stringsAsFactors = F)
  
  
  colnames(t) <- c("nb", "size", "base")
  
  table <- cbind.data.frame("size" = 18:30, 
                            "A" = numeric(13), 
                            "C" = numeric(13), 
                            "G" = numeric(13), 
                            "T" = numeric(13), 
                            "N" = numeric(13))
  
  for(base in unique(t$base)){
    
    tmp <- t[t$base == as.character(base),]
    
    for(i in 1:length(table$size)){
      
      size = table$size[i]
      
      if(nrow(tmp[tmp$size == size,]) > 0){
        
        table[i, base] <- tmp[tmp$size == size, ]$nb
        
      } else {
        
        table[table$size == i, base] <- 0
        
      }
    }
    
  }
  list_table_size_distribution[[index_]] <- table
  write.csv(x = table, file = paste0("table_size_distribution_and_1st_base_",index_,".csv"), 
            quote = F, row.names = F)
}

# for proportion computation

nReads <- unlist(lapply(list_table_size_distribution, function(l){return(sum(l[,2:5]))}))
names(nReads) <- names(list_table_size_distribution)

max_abs <- unlist(lapply(list_table_size_distribution, function(l){return(max(l[,2:5]))}))

max_prop <- unlist(lapply(1:length(list_table_size_distribution), function(i){return(max(list_table_size_distribution[[i]][,2:5]/nReads[i]))}))
names(max_prop) <- index

# plot

for(index_ in index){
  
  i_name <- description_data[which(index == index_)]
  
  t = list_table_size_distribution[[index_]]
  plot(t$size, (t$A/nReads[index_]), type = "l",
       xlab = "read length", 
       ylab = "Proportion", 
       main = paste0("read's 1st nucleotide (", i_name,")"), 
       col = "olivedrab", lwd = 2, 
       ylim = c(0,30), bty = "n")
  lines(t$size, (t$C/nReads[index_]), col = "cornflowerblue", lwd = 2)
  lines(t$size, (t$G/nReads[index_]), col = "firebrick", lwd = 2)
  lines(t$size, (t$T/nReads[index_]), col = "gold", lwd = 2)
  legend("topright", legend = c("A", "C", "G", "T"), 
         col = c("olivedrab", "cornflowerblue", "firebrick" ,"gold"), 
         lty = 1, lwd = 2, bty = "n",y.intersp = 0.7, cex = 0.7)
  
  dev.off() 
}


### Fig S3 B,C & D : densities solo and antisens, SUT, CUT and XUT (small RNAseq)

# prepare list for boxplot

df <- cbind.data.frame(annot_all, 
                       WT = rowMeans(counts_norm[,grep("WT.*densities", colnames(counts_norm))]),
                       xrn1 = rowMeans(counts_norm[,grep("xrn1_[12].*densities", colnames(counts_norm))]),
                       dcr1 = rowMeans(counts_norm[,grep("^dcr1.*densities", colnames(counts_norm))]),
                       xrn1_dcr1 = rowMeans(counts_norm[,grep("xrn1_dcr1.*densities", colnames(counts_norm))]))
df <- df[df$Type %in% c("SUT", "CUT", "XUT"),]
df$is_anti <- df$ID %in% nc.anti.1nt)

df <- melt(data = df, id.vars = c("Type", "is_anti"), measure.vars = c("WT", "xrn1", "dcr1", "xrn1_dcr1"), 
           variable.name = "strain", value.name = "densities")

df$is_anti <- factor(df$is_anti, levels = c(TRUE,FALSE))

list_exprs <- list()
for(t in c("SUT", "CUT", "XUT")){
  i=1
  list_exprs[[t]] <- list()
  for(strain in c("WT", "xrn1", "dcr1", "xrn1_dcr1")){
    for(anti in c(T,F)){
      tmp <- df[df$Type == t & 
                  as.character(df$is_anti) == anti & 
                  as.character(df$strain) == strain,]
      list_exprs[[t]][[i]] <- log2(tmp[tmp$densities != 0,]$densities)
      i <- i+1
    }
  }
}

# boxplot

for(t in c("SUT", "CUT", "XUT")){
  par(mar = c(5,4,4,2))
  boxplot(list_exprs[[t]], outline = F, frame = F,
          col = c("grey80", "grey40"), 
          lwd = 2, medlwd = 3, staplewex = 0, lty = 1, 
          at = c(1,2,4,5,7,8,10,11), 
          xaxt = "n", ylab = "log2 densities", 
          main = paste0("smallRNA signal (22-23nt) for ",t,"s"))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), 
       labels = c(c("WT", "xrn1", "dcr1", "xrn1 dcr1")), tick = F)
  legend("topright", legend = c("antisens", "solo"), 
         text.col = c("grey60", "grey20"), text.width = 3, bty = "n")
}


### Fig 5.A : size antisens CUT in castellii and cerevisiae 

# size asCUT (overlap ORF)

size_list <- list()

tmp <- CUT_antisens_mRNA[!duplicated(CUT_antisens_mRNA$V9),]
size_list[["castellii"]] <- tmp$V5 - tmp$V4 + 1

tmp <- antisens_CUT_ORF_cer[!duplicated(antisens_CUT_ORF_cer$V18),]
size_list[["cerevisiae"]] <- tmp$V14 - tmp$V13 + 1

boxplot(list(size_list$castellii,
             size_list$cerevisiae[size_list$cerevisiae >= 200]),
        outline = F, names =  names(size_list), lwd  = 2, notch = T, medlwd = 3,
        staplelty = 0, ylim = c(0,2500), lty = 1, frame = F, 
        col = rep(rgb(0,176,80, maxColorValue = 255), 2),
        ylab = "Length (nt)", main = "Size of antisens CUT (of ORF, > 200nt)")


### Fig 5.B : size of antisens XUT in castellii and cerevisiae 

# size asXUT (overlap ORF)
size_list <- list()

tmp <- XUT_antisens_mRNA[!duplicated(XUT_antisens_mRNA$V9),]
size_list[["castellii"]] <- tmp$V5 - tmp$V4 + 1

tmp <- antisens_XUT_ORF_cer[!duplicated(antisens_XUT_ORF_cer$V18),]
size_list[["cerevisiae"]] <- tmp$V14 - tmp$V13 + 1

boxplot(size_list, outline = F, col = rep(rgb(255,0,0, maxColorValue = 255), 2),
        staplelty = 0, ylim = c(0,2500), lty = 1, frame = F, lwd = 2, notch = T, medlwd = 3,
        ylab = "Length (nt)", main = "Size of antisens XUT (of ORF)")


### Fig 5.C : size over antisens CUT/ORF in castellii and cerevisiae

overlap_list <- list()
overlap_list[["castellii"]] <- CUT_antisens_mRNA

tmp <- antisens_CUT_ORF_cer
res <- cbind.data.frame(CUT = character(0), gene = character(0), length_orf = numeric(0), length_overlap = numeric(0))
for(cut in unique(tmp$V18)){
  df <- tmp[tmp$V18 == cut,]
  genes_ <- unique(df$V9)
  for(g in genes_){
    df_ <- df[df$V9 == g,]
    l_orf = sum(apply(df_,1,function(v){return(as.numeric(v[5])-as.numeric(v[4])+1)}))
    l_over = sum(df_$V19)
    res <- rbind.data.frame(res, 
                            cbind.data.frame(CUT= cut, gene = g, length_orf = l_orf, length_overlap = l_over))
  }
}
overlap_list[["cerevisiae"]] <- res

size_cut <- annot_all_cer[annot_all_cer$Type == "CUT",]
tmp <- size_cut$Stop - size_cut$Start + 1
names(tmp) <- size_cut$ID
size_cut <- tmp
size_cut <- size_cut[as.character(overlap_list$cerevisiae$CUT)]

boxplot(list(overlap_list$castellii$V19,
             overlap_list$cerevisiae$length_overlap[size_cut >= 200]),
        outline = F, ylim = c(0,2000), col = rep(rgb(0,176,80, maxColorValue = 255), 2),
        names = names(overlap_list), lwd = 2, notch = T,
        staplelty = 0, lty = 1, frame = F, medlwd = 3,
        ylab = "Length (nt)", main = "Size of overlap between ORF and antisens CUT (>200nt)")


### Fig 5.D : size over antisens XUT/ORF in castellii and cerevisiae

overlap_list <- list()

overlap_list[["castellii"]] <- XUT_antisens_mRNA

tmp <- antisens_XUT_ORF_cer
res <- cbind.data.frame(XUT = character(0), gene = character(0), length_orf = numeric(0), length_overlap = numeric(0))
for(xut in unique(tmp$V18)){
  df <- tmp[tmp$V18 == xut,]
  genes_ <- unique(df$V9)
  for(g in genes_){
    df_ <- df[df$V9 == g,]
    l_orf = sum(apply(df_,1,function(v){return(as.numeric(v[5])-as.numeric(v[4])+1)}))
    l_over = sum(df_$V19)
    res <- rbind.data.frame(res, 
                            cbind.data.frame(XUT= xut, gene = g, length_orf = l_orf, length_overlap = l_over))
  }
}
overlap_list[["cerevisiae"]] <- res

boxplot(list(overlap_list$castellii$V19,
             overlap_list$cerevisiae$length_overlap),
        names = names(overlap_list)[1:2], medlwd = 3,
        outline = F, ylim = c(0,2000), col = rep(rgb(255,0,0, maxColorValue = 255), 2),
        staplelty = 0, lty = 1, frame = F, lwd = 2, notch = T,
        ylab = "Length (nt)", main = "Size of overlap between ORF and antisens XUT")
dev.off()

### Fig 5.E : Cumulative coverage of the coding regions by asCUTs and asXUTs in N. castellii and S. cerevisiae 

# N. castellii

# nb nt coding regions

size_coding <- sum(annot_all[annot_all$Type == "gene",]$Stop - annot_all[annot_all$Type == "gene",]$Start + 1)

# nb nt overlap CUT/ coding regions

size_overlap_CUT_asmRNA <- sum(CUT_antisens_mRNA$V19)

# nb nt overlap XUT/ coding regions

size_overlap_XUT_asmRNA <- sum(XUT_antisens_mRNA$V19)

df_cas <- c(size_overlap_CUT_asmRNA/size_coding*100, 
            size_overlap_XUT_asmRNA/size_coding*100)

# S. cerevisiae

# nb nt coding regions

size_coding_cer <- sum(annot_all_cer[annot_all_cer$Type == "gene",]$Stop - annot_all_cer[annot_all_cer$Type == "gene",]$Start + 1)

# nb nt overlap CUT/ coding regions

size_overlap_CUT_asmRNA_cer <- sum(CUT_antisens_mRNA_cer$V19)

# nb nt overlap XUT/ coding regions

size_overlap_XUT_asmRNA_cer <- sum(XUT_antisens_mRNA_cer$V19)

df_cer <- c(size_overlap_CUT_asmRNA_cer/size_coding_cer*100, 
            size_overlap_XUT_asmRNA_cer/size_coding_cer*100)

df <- cbind(df_cas, df_cer) 

# barplot 

barplot(df, beside = T, col = c("grey60", "black"),
        ylab = "% coding region", 
        names.arg = c("asCUT", "asXUT"), 
        main = "coverage coding regions by aslncRNAs",
        ylim = c(0,15))
legend("topright", legend = c("N. castellii", "S. cerevisiae"), fill = c("grey60", "black"), 
       bty = "n", xpd = T)


### Fig S5 : Box-plot for size distribution for the antisens SUTs, CUTs, XUTs identified in this work and for the 170 aslncRNAs previously annotated in N. castellii (Alcid and Tsukiyama 2016)

# prepare data frame 

df <- annot_all[annot_all$Type %in% c("SUT", "CUT", "XUT", "aslncRNA"),]
df$is_anti <- ifelse(df$ID %in% c(SUT_antisens_mRNA$V9, XUT_antisens_mRNA$V9, CUT_antisens_mRNA$V9),
                     'anti', "solo")
df$length <- df$Stop - df$Start + 1
df$Type <- factor(df$Type, levels = c("SUT", "CUT", "XUT", "aslncRNA"))

# boxplot 

boxplot(list(SUT = df[df$is_anti == "anti" & df$Type == "SUT",]$length,
             CUT = df[df$is_anti == "anti" & df$Type == "CUT",]$length,
             XUT = df[df$is_anti == "anti" & df$Type == "XUT",]$length,
             aslncRNA = df[df$Type == "aslncRNA",]$length), 
        notch = T, outline = F, ylim = c(0, 2500), lwd = 2, 
        ylab = "length (nt)", frame = F)
abline(h = 200, col = "red", lty = 2)

