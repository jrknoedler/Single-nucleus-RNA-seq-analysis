output <- "Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_CCAtestallstandardscale_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)


PrimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Primed.rds")
PrimedPOA <- subset(PrimedPOA, idents=c("0","16"), invert=TRUE)
PrimedPOA
MalePOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Intact.rds")
MalePOA <- subset(MalePOA, idents=c("0","15","21","24","28","30","37","35"), invert=TRUE)
MalePOA
UnprimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Prepaperrecluster_2_Unprimed.rds")
UnprimedPOA <- subset(UnprimedPOA, idents=c("0","33"), invert=TRUE)
MalePOA$sex <- "Male"
MalePOA$Hormone <- "Intact"
PrimedPOA$sex <- "Female"
PrimedPOA$Hormone <- "Primed"
UnprimedPOA$sex <- "Female"
UnprimedPOA$Hormone <- "Unprimed"

list <- c(MalePOA,UnprimedPOA, PrimedPOA)
genes.male <- (x=rownames(x=MalePOA))
genes.unprimed <- (x=rownames(x=UnprimedPOA))
genes.primed <- (x=rownames(x=PrimedPOA))

allgenes <- intersect(intersect(genes.male,genes.primed),genes.unprimed)
DEGs <- read.csv("topGO/Total_ByRegion/POA_Genesonly.txt", header=FALSE)
DEGs <- unlist(DEGs)
DEGs <- unique(DEGs)
DEGs.filtered <- intersect(DEGs, allgenes)
genelist <- read.table("ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

DefaultAssay(PrimedPOA) <- "RNA"
DefaultAssay(UnprimedPOA) <- "RNA"
DefaultAssay(MalePOA) <- "RNA"

PrimedPOA <- NormalizeData(PrimedPOA)
UnprimedPOA <- NormalizeData(UnprimedPOA)
MalePOA <- NormalizeData(MalePOA)

PrimedPOA <- FindVariableFeatures(PrimedPOA, selection.method="vst", nfeatures=2000)

UnprimedPOA <- FindVariableFeatures(UnprimedPOA, selection.method="vst", nfeatures=2000)

MalePOA <- FindVariableFeatures(MalePOA, selection.method="vst", nfeatures=2000)

features <- SelectIntegrationFeatures(object.list=list)

features.final <- union(features, DEGs.filtered)
features.final < unlist(features.final)
POA.anchors <- FindIntegrationAnchors(object.list=list, anchor.features = features)

POA.comb <- IntegrateData(anchorset=POA.anchors)

DefaultAssay(POA.comb) <- "integrated"

POA.comb <- ScaleData(POA.comb)
POA.comb <- RunPCA(POA.comb)
POA.comb <- RunUMAP(POA.comb, reduction="pca", dims=1:30)
POA.comb <- FindNeighbors(POA.comb, reduction="pca", dims=1:30)
POA.comb <- FindClusters(POA.comb, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
DefaultAssay(POA.comb) <- "RNA"
POA.comb <- NormalizeData(POA.comb)
pdf(file=paste0(output,"MarkerVln.pdf"), width=32, height=20)
VlnPlot(POA.comb, features=c("Esr1","Pgr","Slc17a6","Gad1","Kiss1"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output,"_Kiss1vln.pdf"), width=20)
VlnPlot(POA.comb, features=c("Kiss1"), split.by="sex", pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(POA.comb, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(POA.comb, file=paste0(output,"CCAscaledata_DEGsonly.rds"))