output <- "Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_Test_CCAtestallstandardscale_DEGenforce"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)


PrimedMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered2_Primed_30pcs_Res1.5_Intact.rds")
MaleMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered3_Intact_30pcs_1.5res_Intact.rds")
UnprimedMeA <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_filtered1v3_30pcs_lowfdr_finalredo.rds")
MaleMeA$sex <- "Male"
MaleMeA$Hormone <- "Intact"
PrimedMeA$sex <- "Female"
PrimedMeA$Hormone <- "Primed"
UnprimedMeA$sex <- "Female"
UnprimedMeA$Hormone <- "Unprimed"

list <- c(MaleMeA,UnprimedMeA, PrimedMeA)
genes.male <- (x=rownames(x=MaleMeA))
genes.unprimed <- (x=rownames(x=UnprimedMeA))
genes.primed <- (x=rownames(x=PrimedMeA))

allgenes <- intersect(intersect(genes.male,genes.primed),genes.unprimed)
DEGs <- read.csv("topGO/Total_ByRegion/MeA_Genesonly.txt", header=FALSE)
DEGs <- unlist(DEGs)
DEGs <- unique(DEGs)
DEGs.filtered <- intersect(DEGs, allgenes)
genelist <- read.table("ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

DefaultAssay(PrimedMeA) <- "RNA"
DefaultAssay(UnprimedMeA) <- "RNA"
DefaultAssay(MaleMeA) <- "RNA"

PrimedMeA <- NormalizeData(PrimedMeA)
UnprimedMeA <- NormalizeData(UnprimedMeA)
MaleMeA <- NormalizeData(MaleMeA)

PrimedMeA <- FindVariableFeatures(PrimedMeA, selection.method="vst", nfeatures=2000)

UnprimedMeA <- FindVariableFeatures(UnprimedMeA, selection.method="vst", nfeatures=2000)

MaleMeA <- FindVariableFeatures(MaleMeA, selection.method="vst", nfeatures=2000)


features <- SelectIntegrationFeatures(object.list=list)

features.final <- union(features, DEGs.filtered)
features.final <- unlist(features.final)
features.final <- setdiff(features.final, genelist)
MeA.anchors <- FindIntegrationAnchors(object.list=list, anchor.features = features.final)

MeA.comb <- IntegrateData(anchorset=MeA.anchors)

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