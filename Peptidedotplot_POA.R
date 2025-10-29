#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)

output <- "Seurat/POA_IndependentAnalysis/POA_peptides_"

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
typemarkers <- read.csv("mouseneuropeptides.csv", header=FALSE)
typemarkers <- unlist(typemarkers)
#mylevels <- c(1,2,3,5,6,10,11,14,15,16,18,19,20,23,24,25,26,0,4,7,8,9,12,13,17,21,22)
#mylevels <- rev(mylevels)
Pct <- DotPlot(mySeurat, features=typemarkers)
data <- Pct$data




pdf(file=paste0(output,"_filtereddotplot_scaledatabluered.pdf"), width=50)
data  %>% ggplot(aes(x=features.plot, y=id, color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"))+ scale_size(limits=c(0,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()