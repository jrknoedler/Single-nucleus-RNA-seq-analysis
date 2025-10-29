#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_Junetest_Superfiltres1.2"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered2_Primed_30pcs_Res1.5_Intact.rds")
MaleMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered3_Intact_30pcs_1.5res_Intact.rds")
UnprimedMeA <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_.0001FDR_1000iter.rds")
UnprimedMeA <- subset(UnprimedMeA, idents=c("3"), invert=TRUE)
#PrimedMeA <- SCTransform(PrimedMeA)
#MaleMeA <- SCTransform(MaleMeA)
#UnprimedMeA <- SCTransform(UnprimedMeA)
MaleMeA$sex <- "Male"
MaleMeA$Hormone <- "Intact"
PrimedMeA$sex <- "Female"
PrimedMeA$Hormone <- "Primed"
UnprimedMeA$sex <- "Female"
UnprimedMeA$Hormone <- "Unprimed"

mySeurat <- merge(MaleMeA, y=c(PrimedMeA,UnprimedMeA), add.cell.ids=c("MaleMeA","PrimedMeA","UnprimedMeA"), project="Merged")
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
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMIvln.pdf"), width=12)
VlnPlot(mySeurat,features=c("nCount_RNA"), pt.size=1)
dev.off() 
pdf(paste0(output,"markervln.pdf"),width=35, height=12)
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Gad1","Cyp19a1","Slc29a4"), pt.size=0, ncol=1)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
props <- prop.table(table(Idents(mySeurat), mySeurat$Hormone), margin=2)
write.csv(props, file=paste0(output, "propperclust.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()