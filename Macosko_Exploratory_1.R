#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNSTMacosko/Macosko_BNST_Exploratory_res1"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

Male.data <- Read10X(data.dir = "Macosko/Male_BNST", gene.column=1)
#Male_Metadata <- read.table("/scratch/users/knoedler/Macosko/Male_Metadata.txt", sep="\t", header=FALSE, row.names=1)
Male <- CreateSeuratObject(counts = Male.data, project = "MaleBNST", min.cells=3, min.features=200)
#Male <- AddMetaData(Male, metadata=Male_Metadata, col.name="Replicate")
Female.data <- Read10X(data.dir="Macosko/Female_BNST", gene.column=1)
Female <- CreateSeuratObject(counts = Female.data, project = "FemaleBNST", min.cells=3, min.features=200)
#Female_Metadata <- read.table("/scratch/users/knoedler/Macosko/Female_Metadata.txt", sep="\t", header=FALSE,row.names=1)
#Female <- AddMetaData(Female, metadata=Female_Metadata, col.name="Replicate")
Male$Hormone <- "Intact"
Male$sex <- "Male"
Female$Hormone <- "Intact"
Female$sex <- "Female"
mySeurat <- merge(Male, y=c(Female), project="MacoskoBNSTMerged")
mySeurat
All_Metadata <- read.table("/scratch/users/knoedler/Macosko/Combined_Metadata.txt", sep="\t", header=FALSE,row.names=1)
mySeurat <- AddMetaData(mySeurat, metadata=All_Metadata, col.name="Replicate")
head(mySeurat[[]])
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
pdf(paste0(output,"_qcVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("nCount_RNA","nFeature_RNA"), ncol=1, pt.size=0, split.by="Replicate")
dev.off()
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("Gad1","Slc17a6","Esr1","Ar","Pgr","Cyp19a1","Tac1"), ncol=1, pt.size=0)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("Gad1","Slc17a6","Esr1","Cyp19a1","Tac1"), ncol=1, pt.size=0)
dev.off()
#pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
#DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cyp19a1","Tac1"))
#dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
#mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
#write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=paste0(output, "_annotated.rds"))

