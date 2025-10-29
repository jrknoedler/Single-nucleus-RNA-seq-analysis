output <- "Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_10pcs_res1_Primed_plots_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Intact <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_10pcs_res1_lowthresh_Primed.rds")
DefaultAssay <- "RNA"
Intact <- NormalizeData(Intact)
pdf(file=paste0(output,"_Esr1Vln.pdf"), width=14, height=14)
VlnPlot(Intact, features=c("Esr1", "Slc17a6","Gad1","Cyp19a1"), pt.size=0, ncol=1)
dev.off()
pdf(file=paste0(output,"_Vln_Markers1.pdf"), width=14, height=40)
VlnPlot(Intact, features=c("Gypa","Prox1","Sytl4","Samd3","Zeb2","Tac2","Slc17a8","Haus4","Cartpt","Mob3b","Cdh23","Lhx5","Ramp1","Satb2","Prdm16","Itga4","Ano1","Avp","Car12","Prkcd"), pt.size=0, ncol=1)
dev.off()
pdf(file=paste0(output,"_Tac1_FeaturePlot.pdf"))
FeaturePlot(Intact, features=c("Tac1"), cols=c("light blue", "dark red"))
dev.off()
