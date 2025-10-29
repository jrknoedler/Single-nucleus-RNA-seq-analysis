#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Pregnancy_Figures/snRNAseqQC_Regressed"
regoutput <- "RegulonsCompiled/POA/Sexmarkersperclust/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)





mySeurat <- readRDS("Pregnancy_Figures/All_snRNAseq_regressed.rds")


pdf(paste0(output,"_UMAP.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()

pdf(paste0(output,"_UMAP_sexlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", shuffle=TRUE, group.by="sex", label=FALSE, cols=c("pink","light blue"))
dev.off()
pc <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("royalblue3", "darkorange2","palevioletred1","chartreuse3"), group.by="Hormone", shuffle=TRUE)
pce <- pc + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())+theme(legend.position="none")


#pcblack <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("yellow", "chartreuse3","deeppink1"), group.by="Hormone", shuffle=TRUE)
#pceblack <- pc + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())+theme(legend.position="none")
#png(paste0(output,"_UMAP_hormonelabel.png"), width=5000,height=5000, units="px")
#pceblack + theme(plot.background=element_rect(fill="black"))
#dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"), width=15,height=15)
pce
dev.off()

p2 <- DimPlot(mySeurat, reduction="umap", shuffle=TRUE, group.by="Region", label=FALSE, cols=c("turquoise2","mediumvioletred","springgreen2","red4"))
p2for <- p2 + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())+theme(legend.position="none")

pdf(paste0(output,"_UMAP_regionlabel.pdf"), width=15,height=15)
p2for
dev.off()
#pdf(paste0(output,"_UMAP_regionlabel_black.pdf"), width=15,height=15)
#p2for + theme(plot.background=element_rect(fill="black"))
#dev.off()

#png(paste0(output,"_UMAP_regionlabel_black.png"), width=5000,height=5000, units="px")
#p2for + theme(plot.background=element_rect(fill="black"))
#dev.off()

p3 <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("yellow4","purple4"), shuffle=TRUE, group.by="Pub")
p3for <- p3+ theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())+theme(legend.position="none")

pdf(paste0(output,"_UMAP_Publabel.pdf"), width=15,height=15)
p3for
dev.off()

#png(paste0(output,"_UMAP_batchlabel_black.png"), width=5000,height=5000, units="px")
#p3for + theme(plot.background=element_rect(fill="black"))
#dev.off()

pdf(paste0(output,"_UMAP_idabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", group.by="orig.ident", label=FALSE)
dev.off()
#saveRDS(mySeurat, file=paste0(output, ".rds"))
Idents(mySeurat) <- "sample"
my_levels <- c("MaleBNST","PrimedBNST","UnprimedBNST","MaleMeA","PrimedMeA","UnprimedMeA","IntactPOA","PrimedPOA","UnprimedPOA","MaleVMH","PrimedVMH","UnprimedVMH")
mySeurat@active.ident <- factor(x=mySeurat@active.ident, levels=my_levels)
pdf(file=paste0(output,"QCUMIplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nCount_RNA"),pt.size=0)
dev.off()
pdf(file=paste0(output,"QCgeneplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nFeature_RNA"), pt.size=0)
dev.off()