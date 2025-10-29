#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_DiffExp/MaleUnprimed_Supergroup"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_filtered.rds")
new.cluster.ids <- c("Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1","Esr1")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
pdf(paste0(output, "_relabeled_UMAP.pdf"))
DimPlot(mySeurat, reduction = "umap", label=TRUE, pt.size=0.5) + NoLegend()
dev.off()
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
Full_Dimorphisms <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", logfc.threshold=0.001, min.pct=0.001, verbose=TRUE)
write.csv(Full_Dimorphisms, file=paste0(output, "_Esr1_default.csv"))
Full_Dimorphisms2 <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", logfc.threshold=0.001, min.pct=0.001, verbose=TRUE, test.use="negbinom")
write.csv(Full_Dimorphisms2, file=paste0(output, "_Esr1_negbinom.csv"))
Full_Dimorphisms3 <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", logfc.threshold=0.001, min.pct=0.001, verbose=TRUE, test.use="poisson")
write.csv(Full_Dimorphisms2, file=paste0(output, "_Esr1_poisson.csv"))
Full_Dimorphisms4 <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", logfc.threshold=0.001, min.pct=0.001, verbose=TRUE, test.use="roc")
write.csv(Full_Dimorphisms2, file=paste0(output, "_Esr1_roc.csv"))
Full_Dimorphisms5 <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", logfc.threshold=0.001, min.pct=0.001, verbose=TRUE, test.use="LR")
write.csv(Full_Dimorphisms2, file=paste0(output, "_Esr1_LR.csv"))
Full_Dimorphisms6 <- FindMarkers(mySeurat, ident.1="Esr1_Male", ident.2="Esr1_Female", min.pct=0, verbose=TRUE, test.use="DESeq2")
write.csv(Full_Dimorphisms2, file=paste0(output, "_Esr1_DESeq2.csv"))
