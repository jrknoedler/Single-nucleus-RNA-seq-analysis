#!/usr/bin/env Rscript

output <- "Seurat/Pregnancy_Full_Dataset_Analysis/Merged_Pregnancy_with_Positive_Regulons_Pseudobulk_Final"

library(Seurat)
#library(pheatmap)
#library(viridis)
#library(RColorBrewer)
#library(dplyr)
#library(tidyverse)
#library(SCENIC)
library(Matrix)
#library(AUCell)
#library(loomR)
library(future)


options(future.globals.maxSize=8000*1024^2)
#path="/scratch/users/knoedler/pySCENIC_Singularity/Subsampled_TFs_OnlyPos_LowNES"
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")
BNST$Region <- "BNST"

ib <- levels(BNST@active.ident)

clustb <- paste0("BNST_",ib)
names(clustb) <- levels(BNST)
BNST <- RenameIdents(BNST, clustb)
BNST$Orig.Clust <- Idents(BNST)

MeA <- readRDS("Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
MeA$Region <- "MeA"

im <- levels(MeA@active.ident)
clustm <- paste0("MeA_",im)
names(clustm) <- levels(MeA)
MeA <- RenameIdents(MeA, clustm)
MeA$Orig.Clust <- Idents(MeA)

POA <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
POA$Region <- "POA"

ip <- levels(POA@active.ident)
clustp <- paste0("POA_",ip)
names(clustp) <- levels(POA)
POA <- RenameIdents(POA, clustp)
POA$Orig.Clust <- Idents(POA)

VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
VMH$Region <- "VMH"


iv <- levels(VMH@active.ident)
clustv <- paste0("VMH_",iv)
names(clustv) <- levels(VMH)
VMH <- RenameIdents(VMH, clustv)
VMH$Orig.Clust <- Idents(VMH)

mySeurat <- merge(x=BNST, y=c(MeA,POA,VMH), add.cell.ids=c("BNST","MeA","POA","VMH"))
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
head(mySeurat[[]])
regulons.master = data.frame()
dim(regulons.master)
#dirs=list.dirs(path, recursive=FALSE)
#dirs

regulons <- read.csv("/scratch/users/knoedler/pySCENIC_Singularity/PseudobulkGRN/GRN_Output/auc_mtx_filtered_2NES_onlyPos_Pseudobulkregulons.csv", check.names=FALSE, header=TRUE, row.names=1)
regulons <- t(regulons)
#regulons.master = rbind(regulons.master, regulons)

#dim(regulons.master)
regs <- as.matrix(regulons)
mySeurat[["Regulons"]] <- CreateAssayObject(data = regs)
mySeurat
VMH <- subset(mySeurat, subset=Region=="VMH")
saveRDS(VMH, file=paste0(output,"VMH_withPosRegulons.rds"))

BNST <- subset(mySeurat, subset=Region=="BNST")
saveRDS(BNST, file=paste0(output,"BNST_withPosRegulons.rds"))

POA <- subset(mySeurat, subset=Region=="POA")
saveRDS(POA, file=paste0(output,"POA_withPosRegulons.rds"))

MeA <- subset(mySeurat, subset=Region=="MeA")
saveRDS(MeA, file=paste0(output,"MeA_withPosRegulons.rds"))


saveRDS(mySeurat, file=paste0(output,"_1stpass_regsonly.rds"))


mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
mySeurat <- FindClusters(mySeurat, resolution=2)
saveRDS(mySeurat, file=paste0(output,"_1stpass_reclustered.rds"))
pdf(paste0(output,"_UMAP.pdf"), width=28, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"), width=42, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"), width=70, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"), width=28, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_regionlabel.pdf"), width=14, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Region")
dev.off()
pdf(paste0(output,"_UMAP_regionsplit.pdf"), width=70, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Region")
dev.off()

mySeurat <- SetIdent(mySeurat, value="Orig.Clust")

pdf(paste0(output,"_Original_Clust_label.pdf"), width=28, height=14)
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
#mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
#write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()



