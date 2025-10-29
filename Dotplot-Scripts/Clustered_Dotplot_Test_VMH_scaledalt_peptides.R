#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_supp7_peptides_25"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
typemarkers <- read.csv("Genelists/Supp6_7_dotplots.csv", header=FALSE)
typemarkers <- unlist(typemarkers)
typemarkers.filtered <- intersect(typemarkers, genes.10x)
head(typemarkers.filtered)
peptidemarkers <- FindAllMarkers(mySeurat, assay="RNA", features=typemarkers.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.223, test.use="MAST")
write.csv(peptidemarkers, file=paste0(output,"Peptidemarkers.csv"))
#mySeurat <- ScaleData(mySeurat)
pdf(file=paste0(output, "Dotplot.pdf"), width=12)
DotPlot(mySeurat, features=typemarkers, dot.min=0.1)
dev.off()
mylevels <- c(1,2,3,5,6,10,11,14,15,16,18,19,20,23,24,25,26,0,4,7,8,9,12,13,17,21,22)
mylevels <- rev(mylevels)
Pct <- DotPlot(mySeurat, features=c(typemarkers))
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
data <- Pct$data




pdf(file=paste0(output,"_filtereddotplot_scaledatabluered.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"))+ scale_size(limits=c(10,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()
