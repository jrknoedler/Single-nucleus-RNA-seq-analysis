output <- "Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_Test_CCAtestallstandardscale_DEGenforce_filtered1"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)



MeA.comb <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_Test_CCAtestallstandardscale_DEGenforceCCAscaledata.rds")
MeA.comb <- subset(MeA.comb, idents=c("19","24","26","30"), invert=TRUE)
DefaultAssay(MeA.comb) <- "integrated"

MeA.comb <- ScaleData(MeA.comb, vars.to.regress="orig.ident")
MeA.comb <- RunPCA(MeA.comb)
MeA.comb <- RunUMAP(MeA.comb, reduction="pca", dims=1:30)
MeA.comb <- FindNeighbors(MeA.comb, reduction="pca", dims=1:30)
MeA.comb <- FindClusters(MeA.comb, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(MeA.comb, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(MeA.comb, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(MeA.comb, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(MeA.comb, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(MeA.comb, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
DefaultAssay(MeA.comb) <- "RNA"
MeA.comb <- NormalizeData(MeA.comb)
pdf(file=paste0(output,"MarkerVln.pdf"), width=32, height=20)
VlnPlot(MeA.comb, features=c("Esr1","Pgr","Slc17a6","Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()
mySeurat.markers <- FindAllMarkers(MeA.comb, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))

saveRDS(MeA.comb, file=paste0(output,"CCAscaledata.rds"))