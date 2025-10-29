output <- "Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_10pcs_res0.8_Unprimed_plots_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Intact <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered2_Unprimed_0.8res_10pcs_res1_Unprimed.rds")
DefaultAssay <- "RNA"
Intact <- NormalizeData(Intact)
pdf(file=paste0(output,"_Vln_Markers1.pdf"), width=14, height=40)
VlnPlot(Intact, features=c("Nfix", "Nxph2","Tmem163","Isl1","Haus4", "Zeb2","Tac2","Cdh23","Lhx5", "Kcnh3","Prdm16","Slc17a8", "Sim1","Kit","Ramp1","Tac1", "Satb2","Eya2","Prkcd","Pax6","Npas1","Crym"), pt.size=0, ncol=1)
dev.off()
pdf(file=paste0(output,"_Tac1_FeaturePlot.pdf"))
FeaturePlot(Intact, features=c("Tac1"), cols=c("light blue", "dark red"))
dev.off()
