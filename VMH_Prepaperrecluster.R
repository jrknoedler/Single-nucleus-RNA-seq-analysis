#!/usr/bin/env Rscript

output <- "Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperrecluster_Sexclude_Malat1regress"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Prepaperrecluster_2_Primed.rds")
PrimedVMH <- subset(PrimedVMH, idents=c("15","16","17"), invert=TRUE)
PrimedVMH
MaleVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Prepaperrecluster_2_Intact.rds")
MaleVMH <- subset(MaleVMH, idents=c("26","0"), invert=TRUE)
MaleVMH
UnprimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Prepaperrecluster_2_Unprimed.rds")
UnprimedVMH <- subset(UnprimedVMH, idents=c("0","1","8","12","14","24"), invert=TRUE)
UnprimedVMH
MaleVMH$sex <- "Male"
MaleVMH$Hormone <- "Intact"
PrimedVMH$sex <- "Female"
PrimedVMH$Hormone <- "Primed"
UnprimedVMH$sex <- "Female"
UnprimedVMH$Hormone <- "Unprimed"
mySeurat <- merge(MaleVMH, y=c(PrimedVMH,UnprimedVMH), add.cell.ids=c("MaleVMH","PrimedVMH","UnprimedVMH"), project="Merged")
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
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Gad1","Cyp19a1","Cckar","Myo1h"),ncol=1, pt.size=0)
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
