#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_CelltypesDraftArcfiltered_regressed"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
a <- AverageExpression(mySeurat, assays=c("RNA"), features=c("Esr1", "Ar","Pgr","Esr2"))
write.csv(a, file=paste0(output,"Esr1perclust.csv"))
typemarkers <- read.table("DotPlotMarkerLists/VMH_Reordered.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)
mylevels <- c(1,2,3,4,6,8,9,10,14,16,17,19,20,21,23,24,25,26,0,5,7,11,12,13,15,18,22,27)
mylevels <- rev(mylevels)
pdf(file=paste0(output,"markerVln.pdf"), width=30, height=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Cckar"), pt.size=0, ncol=1)
dev.off()

pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.25)
dev.off()

Pct <- DotPlot(mySeurat, features=typemarkers)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_rawdotplotscaled.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

# make data square to calculate euclidean distance
mat <- data %>% 
  select(-cluster) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = count) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix

pdf(file=paste0(output,"clusteredDotplot"))
dotplot <- gene_cluster %>% 
  mutate(Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  #filter(pct.exp > 25) %>% 
  ggplot(aes(x=id, y=features.plot, color=avg.exp, size=pct.exp))  + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
dev.off()