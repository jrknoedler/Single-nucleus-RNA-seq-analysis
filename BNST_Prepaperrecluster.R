#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound5_sexclude_malat1exclude_malat1regressres1.2"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Primed.rds")
PrimedBNST <- subset(PrimedBNST, idents=c("2"), invert=TRUE)
PrimedBNST
MaleBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Intact.rds")
MaleBNST <- subset(MaleBNST, idents=c("25","0"), invert=TRUE)
MaleBNST
UnprimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Unprimed.rds")
UnprimedBNST <- subset(UnprimedBNST, idents=c("0"), invert=TRUE)
UnprimedBNST
MaleBNST$sex <- "Male"
MaleBNST$Hormone <- "Intact"
PrimedBNST$sex <- "Female"
PrimedBNST$Hormone <- "Primed"
UnprimedBNST$sex <- "Female"
UnprimedBNST$Hormone <- "Unprimed"
mySeurat <- merge(MaleBNST, y=c(PrimedBNST,UnprimedBNST), add.cell.ids=c("MaleBNST","PrimedBNST","UnprimedBNST"), project="Merged")
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
mySeurat <- FindClusters(mySeurat, resolution=1.2)
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
pdf(paste0(output,"Esr1vln.pdf"),width=35)
VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
dev.off()
pdf(paste0(output,"Tac1vln.pdf"),width=35)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))