#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_DiffExp/MaleUnprimed_NeuronClassdiff_inhib"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_filtered.rds")
new.cluster.ids <- c("Inhibitory_1","Inhibitory_1","Excitatory","Excitatory","Inhibitory_1","Inhibitory_2","Excitatory","Inhibitory_3","Inhibitory_1","Inhibitory_4","Inhibitory_5","Excitatory","Inhibitory_6","Excitatory","Excitatory","Inhibitory_7","Inhibitory_1","Excitatory","Inhibitory_8","Excitatory","Inhibitory_1","Inhibitory_9","Inhibitory_10","Excitatory","Excitatory","Inhibitory_11","Inhibitory_12","Excitatory","Excitatory")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
pdf(paste0(output, "_relabeled_inhib_UMAP.pdf"))
DimPlot(mySeurat, reduction = "umap", label=TRUE, pt.size=0.5) + NoLegend()
dev.off()
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
Excitatory_Dimorphisms <- FindMarkers(mySeurat, ident.1="Inhibitory_1_Male", ident.2="Inhibitory_1_Female", logfc.threshold=0.1, verbose=TRUE)
write.csv(Excitatory_Dimorphisms, file=paste0(output, "_inhibitory_Esr1.csv"))
