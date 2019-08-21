#  *********************************************************
#  *  Script for lncRNA detection from segmentation files  *
#  *********************************************************

# Input : 
#   countFile : RNA expression table 
#   countFileERCC : RNA expression table (for normalization)
#   nReadsWT_1, nReadsWT_2, ... : number of fragment per library, for FPKM computation

# Output : 
#   plot from DESeq (volcano plot, MAplot, etc) & 2D scatter plot were new lncRNA are highlighted
#   table with new XUT, DUT and SUT detected (gff files)
#   Full DESeq results table (fold-change, p-value, etc)
#   CUT have p-value < 0.05, fold-change rrp6/WT > 2 & FPKM in rrp6 strain > 1 FPKM (in all rrp6 strains)
#   SUT have p-value > 0.05 or fold-change mut/WT < 2 & FPKM in WT strain > 1 FPKM (in all WT strains)

# get argumens
args = commandArgs(trailingOnly=TRUE)

countFile = args[1]
countFileCoding = args[2]

nReadWT_1 <- as.numeric(args[3])
nReadWT_2 <- as.numeric(args[4])

nReadRRP6_1 <- as.numeric(args[5])
nReadRRP6_2 <- as.numeric(args[6])

workDIR <- args[7]

nReads <- c(nReadWT_1, nReadWT_2, nReadRRP6_1, nReadRRP6_2)
names(nReads) <- c("WT_1", "WT_2", "rrp6_1", "rrp6_2")

# get counts 
counts <- read.csv(countFile, sep = "\t", comment.char = "#")
colnames(counts) <- c("GeneID", "Chr", "Start", "End", "Strand", "Length", names_samples)

coding <- read.csv(countFileCoding, sep = "\t", comment.char = "#")
colnames(coding) <- c("GeneID", "Chr", "Start", "End", "Strand", "Length", names_samples)

# calculate normalisation coeff (ERCC)
poids = colSums(coding[,names_samples])
coeff = min(poids)/poids

# DESeq
suppressPackageStartupMessages(library(DESeq2))

# prépare tableau de données 
sample_data <- cbind.data.frame(Names = names_samples,
                                Condition = group_samples)

dds <- DESeqDataSetFromMatrix(countData = counts[,names_samples],
                              colData = sample_data,
                              design = ~Condition)
# run DESeq
sizeFactors(dds) <- coeff
dds <- DESeq(dds)

# get raw densities 
# plot WT vs xrn1 and WT vs DICER(raw densities)
for.plot <- counts
# get densities
for.plot.raw <- for.plot[,names_samples]
for.plot.raw <- for.plot.raw/for.plot$Length
# merge samples
for.plot.raw <- cbind.data.frame(WT = log2((for.plot.raw$WT_1 + for.plot.raw$WT_2) / 2), 
                                 rrp6 = log2((for.plot.raw$rrp6_1 + for.plot.raw$rrp6_2) / 2),
                                 stringsAsFactors = F)
for.plot.raw$WT <- replace(for.plot.raw$WT, is.infinite(for.plot.raw$WT), NA)
for.plot.raw$rrp6 <- replace(for.plot.raw$rrp6, is.infinite(for.plot.raw$rrp6), NA)
rownames(for.plot.raw) <- c(as.character(counts$GeneID))

# get normalized densities
for.plot.norm <- for.plot[,names_samples]
for.plot.norm <- data.frame(t(t(for.plot.norm)*coeff))
# densities 
for.plot.norm <- for.plot.norm/for.plot$Length
# merge samples
for.plot.norm <- cbind.data.frame(WT = log2((for.plot.norm$WT_1 + for.plot.norm$WT_2) / 2), 
                                  rrp6 = log2((for.plot.norm$rrp6_1 + for.plot.norm$rrp6_2) / 2))
for.plot.norm$WT <- replace(for.plot.norm$WT, is.infinite(for.plot.norm$WT), NA)
for.plot.norm$rrp6 <- replace(for.plot.norm$rrp6, is.infinite(for.plot.norm$rrp6), NA)

rownames(for.plot.norm) <- c(as.character(counts$GeneID))

# get FPKM
getFPKM <- function(counts_table, gene_length, tot_read){
  res <- matrix(ncol = ncol(counts_table), nrow = nrow(counts_table))
  
  for(i in 1:ncol(counts_table)){
    res[,i] <- (counts_table[,i]/(gene_length*tot_read[i]))*10^9
  }
  
  colnames(res) <- colnames(counts_table)
  
  return(data.frame(res))
}

FPKM <- getFPKM(counts[,names_samples], 
                counts$Length, 
                nReads)

strain ="rrp6"
newType = "CUT"
{
  #plot raw densities
  png(filename = paste0(workDIR, "/2D_plot_WT_vs_",strain,"_raw_densities.png"), width = 540, height = 540)
  par(mar = c(5,4,5,6) + 0.1)
  plot(x = c(min(for.plot.raw$WT, na.rm = T), max(for.plot.raw$WT, na.rm = T)), 
       y = c(min(for.plot.raw[,strain], na.rm = T), max(for.plot.raw[,strain], na.rm = T)), 
       col = "white", 
       xlab = "log2 WT densities", ylab = paste0("log2 ",strain," densities"), 
       main = "raw densities")
  points(for.plot.raw[as.character(counts$GeneID),]$WT, 
         for.plot.raw[as.character(counts$GeneID), strain], 
         col = "black", cex = 1, pch = 21, bg = "grey60",
         xlab = "", ylab = "", 
         main = "")
  abline(0,1)
  legend("topleft", 
         legend = c("new transcripts"), 
         col = "black", pch = 21, pt.bg = "grey60",
         bty = "n", 
         xpd = T)
  dev.off()
  
  #plot norm densities
  png(filename = paste0(workDIR, "/2D_plot_WT_vs_", strain,"_norm_densities.png"), width = 540, height = 540)
  par(mar = c(5,4,5,6) + 0.1)
  plot(x = c(min(for.plot.norm$WT, na.rm = T), max(for.plot.norm$WT, na.rm = T)), 
       y = c(min(for.plot.norm[, strain], na.rm = T), max(for.plot.norm[, strain], na.rm = T)), 
       col = "white", 
       xlab = "log2 WT densities", ylab = paste0("log2 ", strain," densities"), 
       main = "norm densities")
  points(for.plot.norm[as.character(counts$GeneID),]$WT, 
         for.plot.norm[as.character(counts$GeneID), strain], 
         col = "black", cex = 1, pch = 21, bg = "grey60",
         xlab = "", ylab = "", 
         main = "")
  abline(0,1)
  legend("topleft", 
         legend = c("new transcripts"), 
         col = "black", pch = 21, pt.bg = "grey60",
         bty = "n", 
         xpd = T)
  dev.off()
  
  res <- results(dds, alpha = 0.05, contrast = c("Condition", strain, "WT"))
  
  # shrink the log2 fold changes 
  png(paste0(workDIR,"/MA_plot_", strain,"_vs_WT.png"),960,960)
  plotMA(lfcShrink(dds, contrast = c("Condition",strain,"WT"), res = res), ylim = c(-5, 5))
  dev.off()
  
  # volcano plot
  png(paste0(workDIR,"/volcano_plot_", strain, "_vs_WT.png"),960,960)
  plot(res$log2FoldChange, 
       -log10(res$padj), 
       col = ifelse(res$padj < 0.05, "firebrick", "black"),
       bg = ifelse(res$padj < 0.05, "tomato", "grey60"),
       pch = 21, 
       cex = 0.8,
       xlab = "log2 Fold Change", 
       ylab = "-log10 p-value")
  dev.off()
  
  # dispersion estimates
  png(paste0(workDIR,"/dispersion_estimates_", strain, "_vs_WT.png"),960,960)
  plotDispEsts(dds)
  dev.off()
  
  # Define CUT
  
  # fold change > 2 mut/wt and FPKM >= 1 in both mutant samples
  res_zinar <- cbind.data.frame(counts$GeneID, res, stringsAsFactors = F)
  colnames(res_zinar)[1] <- c("ID")
  res_zinar <- replace(res_zinar, is.na(res_zinar), -1) 
  
  tmp <- data.frame(t(t(counts[,names_samples])*coeff))
  foldChange <- ((tmp[,paste0(strain,"_1")] + tmp[,paste0(strain,"_2")]) + 1) / 
    ((tmp[,"WT_1"] + tmp[,"WT_2"]) + 1)
  
  FPKM_gt_1 <- (FPKM[,paste0(strain, "_1")]) >= 1 & (FPKM[,paste0(strain, "_2")]) >=1
  
  new_diff_table <- res_zinar[which(foldChange >= 2 & res_zinar$padj > 0 & res_zinar$padj < 0.05 & FPKM_gt_1),]
  
  tmp <- counts[counts$GeneID %in% new_diff_table$ID, c("GeneID", "Chr", "Start", "End", "Strand")]
  
  if(nrow(tmp) != 0){
    new_diff_gff <- cbind.data.frame(tmp$Chr,
                                     rep("zinar", nrow(tmp)),
                                     rep(newType, nrow(tmp)),
                                     tmp$Start,
                                     tmp$End, 
                                     rep(".", nrow(tmp)),
                                     tmp$Strand, 
                                     rep(".", nrow(tmp)),
                                     paste0("ID=", tmp$GeneID))
    colnames(new_diff_gff) <- c("Chr", "Source", "Type", "Start", "Stop", "", "Strand", "", "ID")
  }
  
  # 2D plot mutant vs WT with new CUT
  png(filename = paste0(workDIR, "/2D_plot_WT_vs_", strain,"_new_", newType,".png"), width = 540, height = 540)
  par(mar = c(5,4,5,6) + 0.1)
  plot(x = c(min(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T), 
             max(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T)), 
       y = c(min(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T), 
             max(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T)), 
       col = "white", 
       xlab = "log2 WT densities", ylab = paste0("log2 ", strain," densities"), 
       main = "norm densities")
  points(for.plot.norm[!(rownames(for.plot.norm) %in% as.character(new_diff_table$ID)),]$WT, 
         for.plot.norm[!(rownames(for.plot.norm) %in% as.character(new_diff_table$ID)), strain], 
         col = "black", bg = "grey60", cex = 1, pch = 21, 
         xlab = "", ylab = "", 
         main = "")
  points(for.plot.norm[as.character(new_diff_table$ID),]$WT, 
         for.plot.norm[as.character(new_diff_table$ID), strain], 
         col = "firebrick", bg = "tomato", cex = 1, pch = 21, 
         xlab = "", ylab = "", 
         main = "")
  abline(0,1)
  legend("topleft", 
         legend = c(paste0("new ", newType)), 
         col = "firebrick", pt.bg = "tomato", pch = 21,
         bty = "n", 
         xpd = T)
  dev.off()
  
  if(nrow(tmp) != 0){
    write.table(new_diff_gff, file = paste0(workDIR,"/new_", newType,".gff"), sep = "\t", col.names = F, row.names = F, quote = F)
  } else {
    file.create(file = paste0(workDIR,"/new_", newType,".gff"))
  }
  write.table(res_zinar, file = paste0(workDIR, "/DESeq_table_full_", strain,"_vs_WT.tab"), sep = "\t", row.names = F, col.names = F, quote = F)
  
}

# Define SUT : new transcripts >= FPKM in WT (both replicates) and not DE rrp vs WT
FPKM_gt_1 <- (FPKM[,"WT_1"]) >= 1 & (FPKM[,"WT_2"]) >=1

source("functions.R")
if(file.info(paste0(workDIR,"/new_CUT.gff"))$size > 0){
  new_CUT <- readGff(paste0(workDIR,"/new_CUT.gff"))
} else { 
  new_CUT <- cbind.data.frame(ID = "")
}

new_SUT <- counts[which(!(counts$GeneID %in% as.character(new_CUT$ID)) & FPKM_gt_1), 
                  c("GeneID", "Chr", "Start", "End", "Strand")]
new_SUT_gff <- cbind.data.frame(new_SUT$Chr,
                                rep("zinar", nrow(new_SUT)),
                                rep("SUT", nrow(new_SUT)),
                                new_SUT$Start,
                                new_SUT$End, 
                                rep(".", nrow(new_SUT)),
                                new_SUT$Strand, 
                                rep(".", nrow(new_SUT)),
                                paste0("ID=", new_SUT$GeneID))
colnames(new_SUT_gff) <- c("Chr", "Source", "Type", "Start", "Stop", "", "Strand", "", "ID")

write.table(new_SUT_gff, file = paste0(workDIR,"/new_SUT.gff"), sep = "\t", col.names = F, row.names = F, quote = F)

{
  # 2D plot mutant vs WT with new SUT
  png(filename = paste0(workDIR, "/2D_plot_WT_vs_", strain,"_new_SUT.png"), width = 540, height = 540)
  par(mar = c(5,4,5,6) + 0.1)
  for.plot.norm <- for.plot.norm[1:nrow(counts),]
  plot(x = c(min(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T), 
             max(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T)), 
       y = c(min(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T), 
             max(for.plot.norm$WT, for.plot.norm[, strain], na.rm = T)), 
       col = "white", 
       xlab = "log2 WT densities", ylab = paste0("log2 ", strain, "  dicer densities"), 
       main = "norm densities")
  points(for.plot.norm[!(rownames(for.plot.norm) %in% gsub("ID=", "", new_SUT_gff$ID)),]$WT, 
         for.plot.norm[!(rownames(for.plot.norm) %in% gsub("ID=", "", new_SUT_gff$ID)), strain], 
         col = "black", bg = "grey60", cex = 1, pch = 21, 
         xlab = "", ylab = "", 
         main = "")
  points(for.plot.norm[gsub("ID=", "", new_SUT_gff$ID),]$WT, 
         for.plot.norm[gsub("ID=", "", new_SUT_gff$ID), strain], 
         col = "firebrick", bg = "tomato", cex = 1, pch = 21, 
         xlab = "", ylab = "", 
         main = "")
  abline(0,1)
  legend("topleft", 
         legend = c("new SUT"), 
         col = "firebrick", pt.bg = "tomato", pch = 21,
         bty = "n", 
         xpd = T)
  dev.off()
  
}
