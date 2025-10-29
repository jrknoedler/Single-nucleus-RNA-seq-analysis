#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound1_sexclude_malat1exclude_malat1regress_filtered4"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
#library(DoubletFinder)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound1_sexclude_malat1exclude_malat1regress_filtered3.rds")
#DefaultAssay(mySeurat) <- "integrated"
mySeurat <- subset(mySeurat, idents=c("32"), invert=TRUE)
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
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
#sweep.res.list <- paramSweep_v3(mySeurat, PCs=1:30, sct=TRUE)
#sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
#bcmvn <- find.pK(sweep.stats)
#bcmvn
#nExp_poi <- round(0.15*nrow(mySeurat@meta.data))
#mySeurat <- doubletFinder_v3(mySeurat, PCs=1:30, pN=0.25, pK=0.14, nExp=nExp_poi, reuse.pANN=FALSE, sct=TRUE)
#head(mySeurat[[]])
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
pdf(paste0(output,"Esr1vln.pdf"),width=35)
VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
dev.off()
pdf(paste0(output,"Tac1vln.pdf"),width=35)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
dev.off()
pdf(paste0(output,"markervln.pdf"),width=35, height=12)
VlnPlot(mySeurat, features=c("Slc17a6","Gad1","Cyp19a1"), pt.size=0, ncol=1)
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
#pdf(file=paste0(output,"doublettest.pdf"))
#DimPlot(mySeurat, group.by="DF.classifications_0.25_0.14_1476")
#dev.off()


