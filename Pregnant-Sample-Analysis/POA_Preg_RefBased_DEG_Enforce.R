#!/usr/bin/env Rscript
rm(list=ls())

library(scclusteval)
library(MAST)
library(future)
library(Matrix)

library(tidyverse)
library(patchwork)
library(Seurat)
library(dplyr)


options(future.globals.maxSize=8000*1024^2)

output <- "/scratch/users/knoedler/Seurat/POA_Preg_RefBased"

##Load dataset 1
Published <- readRDS("/scratch/groups/nirao/shared1/POA_Published.rds")
Published@meta.data$NewIntegration <- "Published"


##Load dataset 2
Pregnant <- readRDS("/scratch/groups/nirao/shared1/PregPOA_1stPass.rds")


Pregnant <- RenameCells(object = Pregnant, add.cell.id = "Pregnant")
Pregnant@meta.data$sex <- "Female"
Pregnant@meta.data$Hormone <- "Pregnant"
Pregnant@meta.data$NewIntegration <- "Pregnant"

#Set default assay to SCT
DefaultAssay(Published) <- "SCT"
DefaultAssay(Pregnant) <- "SCT"

##Create list of objects to be integrated
list <- c(Published, Pregnant)

##Get all genes present in Published data
genes.published <- (x=rownames(x=Published@assays$SCT))

##Get all genes present in Pregnant data
genes.pregnant <- (x=rownames(x=Pregnant@assays$SCT))

##Get only genes present in BOTH datasets
allgenes <- intersect(genes.published, genes.pregnant)

##Get relevant list of DEGs from TRAPseq experiment (e.g., all genes that are DEGs vs pregnant animals). You can also combine multiple lists here.
##Note that you may need to change the command based on file format (e.g. "read.table" vs "read.cvs"), and check headers etc to make
##sure you get everything.
DEGs1 <- read.csv("/scratch/groups/nirao/shared1/POA_AllvFemale1.5_2.txt", header = FALSE)

##reformat object
DEGs1 <- unlist(DEGs1)


##Remove any duplicates on list (like if you combined four partially overlapping lists)
DEGs <- unique(DEGs1)

##Get all DEGs that are detected in snRNAseq datset. 
DEGs.filtered <- intersect(DEGs, allgenes)

##Get genes to exclude from clustering (Xist, Tsix and y-chromosome genes)
exclude <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE) 
exclude <- unlist(exclude)
exclude <- as.matrix(exclude)

##features 
list <- lapply(X = list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
features.degs <- union(features, DEGs.filtered)

#remove any y-chromosome or Xist/Tsix stuff if Seurat decided to use them as PCA features.
#The purpose of all this filtering is simple: if we tell Seurat to use our list of features instead
#of what it thinks is best, every gene we use needs to actually have been detected (otherwise it'll
#throw an error and crash)
features.final <- setdiff(features.degs, exclude)
features.final < unlist(features.final)

list <- PrepSCTIntegration(list, anchor.features = features.final, verbose=TRUE)
list <- lapply(X = list, FUN = RunPCA, features = features.final)

anchors <- FindIntegrationAnchors(object.list = list, reference = c(1),
                                  normalization.method = "SCT", anchor.features = features.final,
                                  reduction = "rpca", dims = 1:50, verbose = TRUE)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT",
                            dims = 1:50, verbose = FALSE)

integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction ="pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution=1.2)

pdf(paste0(output,"_UMAP.pdf"))
DimPlot(integrated, reduction="umap", label=TRUE)
dev.off()

pdf(paste0(output,"_UMAP_sexsplit.pdf"), width = 8.5, height = 4)
DimPlot(integrated, reduction="umap", label=TRUE, split.by="sex")
dev.off()

pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(integrated, reduction="umap", label=TRUE, group.by="sex")
dev.off()

pdf(paste0(output,"_UMAP_hormonesplit.pdf"), width = 15, height = 4)
DimPlot(integrated, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()

pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(integrated, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()

##Find marker genes 
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)

pdf(file=paste0(output,"_MarkerVln.pdf"), width=32, height=20)
VlnPlot(integrated, features=c("Esr1","Pgr","Slc17a6","Gad1"), pt.size=0, ncol=1)
dev.off()


mySeurat.markers <- FindAllMarkers(integrated, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))


saveRDS(integrated, file=paste0(output,".rds"))


##Jaccard Index
Published@meta.data$PublishedCluster <- Idents(Published)
Integration <- integrated
Integration@meta.data$IntegrationCluster <- Idents(Integration)


##Reference
##https://crazyhottommy.github.io/EvaluateSingleCellClustering/mixture_tidy_idents.html
##https://github.com/crazyhottommy/scclusteval/blob/master/R/scclusterboot.R


##Jaccard Similarity Heatmap
pdf(file=paste0(output,"_PairWiseJaccardSetsHeatmap.pdf"), width=27, height=25)
PairWiseJaccardSetsHeatmap(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                           set_names(Integration@meta.data$IntegrationCluster, nm=colnames(Integration)),
                           show_row_dend = F, show_column_dend = F, cluster_row = F, cluster_column =F,
                           col_low = "#FFFFFF", col_high = "#004D40")
dev.off()


##Jaccard Index Calculation
PairWiseJaccardSets <- PairWiseJaccardSets(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                                           set_names(Integration@meta.data$IntegrationCluster, nm=colnames(Integration)))
write.csv(PairWiseJaccardSets, file = paste0(output, "_PairWiseJaccardSets.csv"))








