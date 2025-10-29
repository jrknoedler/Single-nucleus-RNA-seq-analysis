#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper1stpass/20pcs/BNST_IntactvPrimed_Cluster_DESeq2"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_20pcs.rds")

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)

datar <- round(data)

RoundedSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone

RoundedSeurat <- AddMetaData(RoundedSeurat, metadata=idents, col.name="Ident")
RoundedSeurat <- AddMetaData(RoundedSeurat, metadata=hormone, col.name="Hormone")
Idents(RoundedSeurat) <- "Ident"
head(RoundedSeurat[[]])


RoundedSeurat$celltype.Hormone <- paste(Idents(RoundedSeurat), RoundedSeurat$Hormone, sep="_")
RoundedSeurat$celltype <- Idents(RoundedSeurat)
Idents(RoundedSeurat) <- "celltype.Hormone"
resultsFC <- cbind(0)
resultsFC <- data.frame(resultsFC)
resultsP <- cbind(0)
resultsP <- data.frame(resultsP)
for (i in 0:31){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(RoundedSeurat, slot="counts", ident.1 = ident1, ident.2=ident2, verbose=TRUE, test.use="DESeq2")
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
result_i_FC <- sex.dimorphism %>% select(avg_logFC)
result_i_P <- sex.dimorphism %>% select(p_val_adj) 
resultsFC <- merge(resultsFC, result_i_FC, by=row.names)
resultsP <- merge(resultsP, result_i_P, by=row.names)

})
}
write.csv(resultsFC, file=paste0(output,"compiledFCperclust.csv"))
write.csv(resultsP, file=paste0(output, "compiledPperclust.csv"))