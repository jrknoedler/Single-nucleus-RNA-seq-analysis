directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_CelltypesDraft_Malat1regress_Esr1filtered_sexclude_res1.2alt2reordered_reclustered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)
library(patchwork)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
#mySeurat <- NormalizeData(mySeurat)
typemarkers <- read.table("DotPlotMarkerLists/BNST_Reordered_2.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)


mylevels <- c(6,17,22,34,0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33)
mylevels <- rev(mylevels)
Pct <- DotPlot(mySeurat, features=typemarkers, dot.min=0.2)
data <- Pct$data
BNSTp <- data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100), breaks=seq(25,100,5)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())


MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
MeAtypemarkers <- read.table("DotPlotMarkerLists/MeAReorganized_Reordered_2.txt")
MeAtypemarkers <- unlist(MeAtypemarkers)
MeAtypemarkers <- as.matrix(MeAtypemarkers)
MeA <- NormalizeData(MeA)
MeA <- ScaleData(MeA)
PctMeA <- DotPlot(MeA, features=MeAtypemarkers, dot.min=0.2)
MeAdat <- PctMeA$data
MeAlevels <- c(0,3,8,12,13,14,16,17,20,21,22,28,30,31,32,33,1,2,5,6,7,10,11,15,18,19,23,24,26,27,29,34,4,9,25)
MeAlevels <- rev(MeAlevels)
MeAp <- MeAdat %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
POA <- NormalizeData(POA)
POA <- ScaleData(POA)
POAtypemarkers <- read.table("DotPlotMarkerLists/POA_Reordered.txt")
POAtypemarkers <- unlist(POAtypemarkers)
POAtypemarkers <- as.matrix(POAtypemarkers)
POAlevels=c(0,2,7,8,9,15,16,18,22,28,30,36,37,38,1,3,4,5,6,10,11,12,13,14,17,19,20,21,23,24,25,26,27,29,31,32,33,34,35)
POAlevels <- rev(POAlevels)
PctPOA <- DotPlot(POA, features=POAtypemarkers, dot.min=0.2)
POAdat <- PctPOA$data
POAp <- POAdat %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
VMH <- NormalizeData(VMH)
VMHtypemarkers <- read.table("DotPlotMarkerLists/VMH_Reordered_2.txt")
VMHtypemarkers <- unlist(VMHtypemarkers)
VMHtypemarkers <- as.matrix(VMHtypemarkers)
VMH <- ScaleData(VMH)
VMHtypemarkers
VMHlevels <- c(1,2,3,4,6,8,9,10,14,16,17,19,20,21,23,24,25,26,0,5,7,11,12,13,15,18,22,27)
VMHlevels <- rev(VMHlevels)
PctVMH <- DotPlot(VMH, features=VMHtypemarkers, dot.min=0.2)
VMHdat <- PctVMH$data
VMHp <- VMHdat  %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

pdf(file="Seurat/AllDotplots_test8.pdf", width=82, height=8)
BNSTp+MeAp+POAp+VMHp+plot_layout(ncol=4, widths=c(1,1,1.09,0.84))
dev.off()
#pdf(file="Seurat/AllDotplots_test2.pdf")
#BNSTp+MeAp+POAp+VMHp+plot_layout(ncol=1, widths=c(1,1,1.09,0.84))
#dev.off()
#pdf(file="Seurat/AllDotplots_test3.pdf", width=30, height=40)
#BNSTp+MeAp+POAp+VMHp+plot_layout(ncol=1)
#dev.off()