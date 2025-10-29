output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_CCAtestallstandardscale_DEGenforce"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)


PrimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_40pcs_res1.5_lowthresh_Primed.rds")
MaleBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_40pcs_res1_lowthresh_res1.5_Intact.rds")
UnprimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered2_Unprimed_30pcs_nothresh_res1__Unprimed.rds")
MaleBNST$sex <- "Male"
MaleBNST$Hormone <- "Intact"
PrimedBNST$sex <- "Female"
PrimedBNST$Hormone <- "Primed"
UnprimedBNST$sex <- "Female"
UnprimedBNST$Hormone <- "Unprimed"

list <- c(MaleBNST,UnprimedBNST, PrimedBNST)
genes.male <- (x=rownames(x=MaleBNST))
genes.unprimed <- (x=rownames(x=UnprimedBNST))
genes.primed <- (x=rownames(x=PrimedBNST))

allgenes <- intersect(intersect(genes.male,genes.primed),genes.unprimed)
DEGs <- read.csv("topGO/Total_ByRegion/BNST_Genesonly.txt", header=FALSE)
DEGs <- unlist(DEGs)
DEGs <- unique(DEGs)
DEGs.filtered <- intersect(DEGs, allgenes)
genelist <- read.table("ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

DefaultAssay(PrimedBNST) <- "RNA"
DefaultAssay(UnprimedBNST) <- "RNA"
DefaultAssay(MaleBNST) <- "RNA"

PrimedBNST <- NormalizeData(PrimedBNST)
UnprimedBNST <- NormalizeData(UnprimedBNST)
MaleBNST <- NormalizeData(MaleBNST)

PrimedBNST <- FindVariableFeatures(PrimedBNST, selection.method="vst", nfeatures=2000)

UnprimedBNST <- FindVariableFeatures(UnprimedBNST, selection.method="vst", nfeatures=2000)

MaleBNST <- FindVariableFeatures(MaleBNST, selection.method="vst", nfeatures=2000)

features <- SelectIntegrationFeatures(object.list=list)

features.final <- union(features, DEGs.filtered)
features.final <- unlist(features.final)
features.final <- setdiff(features.final, genelist)
BNST.anchors <- FindIntegrationAnchors(object.list=list, anchor.features = features.final)

BNST.comb <- IntegrateData(anchorset=BNST.anchors)

DefaultAssay(BNST.comb) <- "integrated"

BNST.comb <- ScaleData(BNST.comb, vars.to.regress="orig.ident")
BNST.comb <- RunPCA(BNST.comb)
BNST.comb <- RunUMAP(BNST.comb, reduction="pca", dims=1:30)
BNST.comb <- FindNeighbors(BNST.comb, reduction="pca", dims=1:30)
BNST.comb <- FindClusters(BNST.comb, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(BNST.comb, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(BNST.comb, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(BNST.comb, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(BNST.comb, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(BNST.comb, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
DefaultAssay(BNST.comb) <- "RNA"
BNST.comb <- NormalizeData(BNST.comb)
pdf(file=paste0(output,"MarkerVln.pdf"), width=40, height=20)
VlnPlot(BNST.comb, features=c("Esr1","Pgr","Slc17a6","Gad1","Cyp19a1","Tac1"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output,"_Tac1vln.pdf"), width=20)
VlnPlot(BNST.comb, features=c("Tac1"), split.by="sex", pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(BNST.comb, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(BNST.comb, file=paste0(output,"CCAscaledata_DEGsonly.rds"))