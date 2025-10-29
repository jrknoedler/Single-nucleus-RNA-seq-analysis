#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_RegulonAssay"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(loomR)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered2_sexclude_malat1regress_30pcsres1.2FINAL.rds")
regulons <- read.csv("/scratch/users/knoedler/pySCENIC_Singularity/BNST_DefaultNES_Final/auc_mtx_filtered.csv", check.names=FALSE, header=TRUE, row.names=1)
regulons <- t(regulons)
mySeurat[["Regulons"]] <- CreateAssayObject(data = regulons)

head(mySeurat[[]])
head(mySeurat[["Regulons"]])

Degulators <- read.table("RegulonsCompiled/BNST/AllDegulators.txt")
Degulators <- unlist(Degulators)
Degulators <- as.matrix(Degulators)
Degulators
mySeurat <- ScaleData(mySeurat, assay="Regulons", features=Degulators)
AUC <- FindAllMarkers(mySeurat, assay="Regulons", features=Degulators, test.use="roc", logfc.threshold=0, min.pct=0, only.pos=TRUE)
write.csv(AUC, file=paste0(output,"_regulonROCbyclust.csv"))


pdf(paste0(output,"_RegulonROC_heatmap.pdf"))
top10 <- AUC %>% group_by(cluster) %>% top_n(n=2, wt=power)
DoHeatmap(mySeurat, features=top10$gene, raster=FALSE, assay="Regulons") 
dev.off()

DefaultAssay(mySeurat) <- "Regulons"

pdf(file=paste0(output,"_Tcf4UMAP.pdf"), height=60)
FeaturePlot(mySeurat, features=c("Tcf4(+)", "Otp(+)","Arx(+)","Cux1(+)","Nrf1(+)"), slot="scale.data", reduction="umap", cols=c("light blue", "red"))
dev.off()



pdf(file=paste0(output,"_ArxUMAP.pdf"))
FeaturePlot(mySeurat, features=c("Arx(+)"), slot="scale.data", reduction="umap", cols=c("light blue", "red"))
dev.off()


pdf(file=paste0(output,"_OtpUMAP.pdf"))
FeaturePlot(mySeurat, features=c("Otp(+)"), slot="scale.data", reduction="umap", cols=c("light blue", "red"))
dev.off()

mySeurat <- RunPCA(mySeurat, assay="Regulons", features=Degulators, dims=1:30)

mySeurat <- BuildClusterTree(mySeurat, assay="Regulons", dims=1:30, features=Degulators, reorder=TRUE)
pdf(file=paste0(output, "Degulatorregulon_tree.pdf"))
PlotClusterTree(mySeurat)
dev.off()

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

pdf(file=paste0(output, "Arx_VlnTargets.pdf), width=20, height=20)
VlnPlot(mySeurat, features=c("Arx","St18","Nfix","Calb1"), ncol=1, idents=c("14","3","24","17","22"), split.by="Hormone", cols=c("light blue","pink","light green"))
dev.off()

pdf(file=paste0(output, "OtpTargets.pdf"), width=20)
VlnPlot(mySeurat, features=c("Calcrl","Cd36"), idents=c("6","34","23"), ncol=1, split.by="Hormone", cols=c("light blue","pink","light green"))
dev.off()

