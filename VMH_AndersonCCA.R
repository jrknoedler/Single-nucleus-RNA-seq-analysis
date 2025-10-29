#!/usr/bin/env Rscript

output <- "Seurat/VMH_Anderson/Anderson_JRK_Integrated_firstpass"
library(Seurat)
#library(cowplot)
#library(reticulate)
#library(ggplot2)
library(dplyr)
#library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)
snEsr1 <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
macosko <- readRDS("Seurat/VMH_Anderson/Anderson_NonneuronalFilt5.rds")
snEsr1$type <- "Esr1lineage"
macosko$type <- "AllMeA"
snEsr1 <- SCTransform(snEsr1)
DefaultAssay(snEsr1) <- "SCT"
DefaultAssay(macosko) <- "SCT"
genes.sn <- (x=rownames(x=snEsr1))
genes.macosko <- (x=rownames(x=macosko))

allgenes <- intersect(genes.sn,genes.macosko)


list <- c(macosko, snEsr1)
#VMH.features <- read.csv("Seurat/VMH_IndependentAnalysis/VMH_integrationfeatures.txt", header=FALSE)
#VMH.features <- unlist(VMH.features)
#VMH.features <- unique(VMH.features)
#anchors.filtered <- intersect(VMH.features, allgenes)
#genelist <- read.table("ExcludeIDs.txt", header=FALSE)
#genelist <- unlist(genelist)
#genelist <- as.matrix(genelist)
#anchors.final <- setdiff(anchors.filtered, genelist)
Int.features <- SelectIntegrationFeatures(object.list=list, nfeatures=3000)
Sample.list <- PrepSCTIntegration(list, anchor.features=Int.features, verbose=TRUE)
Int.anchors <- FindIntegrationAnchors(Sample.list, normalization.method="SCT", anchor.features=Int.features, verbose=TRUE)
Integrated <- IntegrateData(anchorset=Int.anchors, normalization.method="SCT", verbose=FALSE)
Integrated <- RunPCA(Integrated)
Integrated <- RunUMAP(Integrated, reduction="pca", dims=1:30)
Integrated <- FindNeighbors(Integrated, reduction ="pca", dims=1:30)
Integrated <- FindClusters(Integrated, resolution=1.2)
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
saveRDS(Integrated, file=(paste0(output, ".rds")))