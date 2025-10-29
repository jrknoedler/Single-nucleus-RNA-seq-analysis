#!/usr/bin/env Rscript

output <- "Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperrecluster_Sexclude_Malat1regress"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Primed.rds")
PrimedPOA <- subset(PrimedPOA, idents=c("0","16"), invert=TRUE)
PrimedPOA
MalePOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Intact.rds")
MalePOA <- subset(MalePOA, idents=c("0","15","21","24","28","30","37","35"), invert=TRUE)
MalePOA
UnprimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Unprimed.rds")
UnprimedPOA <- subset(UnprimedPOA, idents=c("0","33"), invert=TRUE)
UnprimedPOA
MalePOA$sex <- "Male"
MalePOA$Hormone <- "Intact"
PrimedPOA$sex <- "Female"
PrimedPOA$Hormone <- "Primed"
UnprimedPOA$sex <- "Female"
UnprimedPOA$Hormone <- "Unprimed"
mySeurat <- merge(MalePOA, y=c(PrimedPOA,UnprimedPOA), add.cell.ids=c("MalePOA","PrimedPOA","UnprimedPOA"), project="Merged")
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")

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
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"Markerln.pdf"),width=35, height=35)
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Gad1","Thbs4","Kiss1","Slc18a2"), ncol=1, pt.size=0)
dev.off()
pdf(paste0(output,"_UMIvln.pdf"), width=30)
VlnPlot(mySeurat,features=c("nCount_RNA"), pt.size=1)
dev.off() 
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))