#!/usr/bin/env Rscript

output <- "Seurat/Dulac_Integrated/Dulac_JRK_Integrated_nonneuronalfiltreanchoredres1.5dim40"
library(Seurat)
#library(cowplot)
#library(reticulate)
#library(ggplot2)
library(dplyr)
#library(MAST)
#library(future)
#options(future.globals.maxSize=1000*1024^2)
genelist <- read.table("ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

snEsr1 <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
snEsr1
hvg1 <- snEsr1@assays$SCT@var.features
hvg1 <- unlist(hvg1)
hvg1 <- as.matrix(hvg1)

macosko <- readRDS("Dulac_Annotated.rds")
macosko
macosko <- subset(macosko, idents=c("Non-Neuronal"), invert=TRUE)
DefaultAssay(snEsr1) <- "SCT"
macosko <- SCTransform(macosko)
hvg2 <- macosko@assays$SCT@var.features
hvg2 <- unlist(hvg2)
hvg2 <- as.matrix(hvg2)
hvgtotal <- union(hvg1, hvg2) 
snEsr1$type <- "Esr1lineage"
macosko$type <- "AllVMH"
DefaultAssay(snEsr1) <- "SCT"
DefaultAssay(macosko) <- "SCT"
genes.sn <- (x=rownames(x=snEsr1))
genes.macosko <- (x=rownames(x=macosko))

allgenes <- intersect(genes.sn,genes.macosko)


list <- c(macosko, snEsr1)
Int.features <- SelectIntegrationFeatures(object.list=list, nfeatures=3000)
features.final <- union(Int.features,hvgtotal)
features.final <- unlist(features.final)
features.final2 <- intersect(features.final, allgenes)
features.final3 <- setdiff(features.final2,genelist)
Sample.list <- PrepSCTIntegration(list, anchor.features=features.final3, verbose=TRUE)
Int.anchors <- FindIntegrationAnchors(Sample.list, normalization.method="SCT", anchor.features=features.final3, verbose=TRUE)
Integrated <- IntegrateData(anchorset=Int.anchors, normalization.method="SCT", verbose=FALSE)
Integrated <- RunPCA(Integrated)
Integrated <- RunUMAP(Integrated, reduction="pca", dims=1:40)
Integrated <- FindNeighbors(Integrated, reduction ="pca", dims=1:40)
Integrated <- FindClusters(Integrated, resolution=1.5)
props <- prop.table(table(Idents(Integrated), Integrated$type), margin = 2)
write.csv(props, file=paste0(output,"Proptable.csv"))
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(Integrated, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_modesplit.pdf"))
DimPlot(Integrated, reduction="umap", label=TRUE, split.by="type")
dev.off()
pdf(paste0(output,"_UMAP_modelabel.pdf"))
DimPlot(Integrated, reduction="umap", label=TRUE, group.by="type")
dev.off()
#DefaultAssay(Integrated) <- "RNA"
#Integrated <- NormalizeData(Integrated)
#mySeurat.markers <- FindAllMarkers(Integrated, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
#write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#pdf(paste0(output,"_Cckarvln.pdf"), width=20)
#VlnPlot(Integrated, features=c("Cckar"), split.by="sex", pt.size=0)
#dev.off()
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
Markers <- read.table("DotPlotMarkerLists/POA_ReorderedFinal_CPM.txt", header=FALSE)
Markers <- unlist(Markers)
pdf(file=paste0(output,"Vlnplottest.pdf"), width=30, height=200)
VlnPlot(Integrated, features=Markers, split.by="type", ncol=1, pt.size=0)
dev.off()
saveRDS(Integrated, file=(paste0(output, ".rds")))