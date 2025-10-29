#!/usr/bin/env Rscript


output <- "Seurat/POA_PregIntegration/Preg_Integration_Testrun1"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

###Load objects to be integrated
#Load reference dataset from Knoedler et al 2022
Ref <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
#Rename Ref clusters based on paper to avoid confusion (numbering starts at 0 by default)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(Ref)
Ref <- RenameIdents(Ref, new.cluster.ids)
#Rerun UMAP using original parameters to get the umap model
Ref <- RunUMAP(Ref, reduction="pca", dim=1:30, return.model=TRUE)

#Load latest filtered pregnant POA single cell data set
Preg <- readRDS("Seurat/PregPOA_Exploratory/PregPOA_Filtered3.rds")
DefaultAssay(Preg) <- "RNA"
Preg <- NormalizeData(Preg)
#Get all genes expressed in each object
Genes.Ref <- (x=rownames(x=Ref))
unlist(Genes.Ref)

Genes.Preg <- (x=rownames(x=Preg))
unlist(Genes.Preg)

#Get all genes present in both datasets
Genes.total <- intersect(Genes.Ref, Genes.Preg)

#Get variable features from each object
Features.Ref <- Ref@assays$SCT@var.features
Features.Preg <- Preg@assays$SCT@var.features

#Get the union of variable features
Features.Total <- union(Features.Ref, Features.Preg)

#Filter genes absent from one data set or another
Features.filtered <- intersect(Features.Total, Genes.total)

#Find integration anchors based on 
Ref.anchors <- FindTransferAnchors(reference=Ref, query=Preg, dims=1:30, reference.reduction="pca", features=Features.filtered)
predictions <- TransferData(anchorset=Ref.anchors, refdata=Ref@active.ident, dims=1:30)
Preg <- AddMetaData(Preg, metadata=predictions)

Preg <- MapQuery(anchorset=Ref.anchors, reference=Ref, query=Preg, refdata=Ref@active.ident, reference.reduction="pca", reduction.model="umap")
saveRDS(Preg, file=paste0(output, "IntegrationTest_refmapped.rds"))
p1 <- DimPlot(Ref, reduction="umap", label=TRUE)

p2 <- DimPlot(Preg, reduction="ref.umap", group.by="predicted.celltype", label=TRUE)
pdf(file=paste0(output, "_InitiaMAP.pdf"), height=7, width=21)
p1+p2
dev.off()



