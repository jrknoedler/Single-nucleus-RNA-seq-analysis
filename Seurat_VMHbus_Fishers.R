#!/usr/bin/env Rscript
output <- "Seurat/VMH_BusMergeTest/Fishers_Test/NoCastrate_filtered1_ExcludeDE_"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_ExcludeDE.rds")
mySeurat$celltype.status <- paste(Idents(mySeurat), mySeurat$status, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.status"
genelist <- read.table("Genelists/VMH_allDEnew.txt")
unlist(genelist)
genelist <- as.matrix(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
for (i in 0:34){
results <- paste0(i,"results")
results =NULL
try({
for (g in filtered.genelist){
try({
identMale <- paste0(i,"_IntactMale")
identFemale <- paste0(i,"_PrimedFemale")
subset1 <- subset(mySeurat, idents=identMale)
subset2 <- subset(mySeurat, idents=identFemale)
g_Group1pos <- sum(GetAssayData(subset1, slot="data")[g,]>0)
g_Group1neg <- sum(GetAssayData(subset1, slot="data")[g,]==0)
g_Group2pos <- sum(GetAssayData(subset2, slot="data")[g,]>0)
g_Group2neg <- sum(GetAssayData(subset2, slot="data")[g,]==0)
g_testor <- rbind(c(g_Group1pos, g_Group1neg), c(g_Group2pos, g_Group2neg))
g_chi2 = chisq.test(g_testor, correct=F)
g_Group1pct <- sum(GetAssayData(object=subset1, slot="data")[g,]>0)/nrow(subset1@meta.data)
g_Group2pct <- sum(GetAssayData(object=subset2, slot="data")[g,]>0)/nrow(subset2@meta.data)
results_g <- data.frame(g, g_Group1pct, g_Group2pct, g_chi2$p.value)
results = rbind(results, results_g)
})
}
})
write.csv(results, file=paste0(output,i,"_fishers_stats.csv"))
}
