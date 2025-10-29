output <- "Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1_15pcs_filtered4_res1_Intact_plots_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Intact <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1_15pcs_filtered4_res1_Intact.rds")
Intact
DefaultAssay <- "RNA"
Intact <- NormalizeData(Intact)
pdf(file=paste0(output,"_Esr1Vln.pdf"), width=14, height=14)
VlnPlot(Intact, features=c("Esr1", "Slc17a6","Gad1","Gad2"), pt.size=0, ncol=1)
dev.off()
pdf(file=paste0(output,"_Tac1_FeaturePlot.pdf"))
FeaturePlot(Intact, features=c("Tac1"), cols=c("light blue", "dark red"))
dev.off()
