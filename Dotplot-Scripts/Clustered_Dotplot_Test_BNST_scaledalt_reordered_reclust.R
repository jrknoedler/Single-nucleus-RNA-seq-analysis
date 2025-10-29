#!/usr/bin/env Rscript

sampleID <- "MalePOA"
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


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis//BNST_Fullmerge_Test_prepaperreclusterRound1_sexclude_malat1exclude_malat1regress_filtered2.rds")
DefaultAssay(mySeurat) <- "RNA"
#mySeurat <- NormalizeData(mySeurat)
typemarkers <- read.table("DotPlotMarkerLists/BNST_Reordered.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.2)
dev.off()

mylevels <- c(5,6,13,0,1,2,3,4,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)
mylevels <- rev(mylevels)
Pct <- DotPlot(mySeurat, features=typemarkers)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_filtereddotplot_scaledata.pdf"), width=30)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point(position=position_dodge(width=0.5)) + scale_colour_gradientn(colours = viridis(500), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank()) + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"_filtereddotplot_scaledatamagma.pdf"), width=30)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point(position=position_dodge(width=1)) + scale_colour_gradientn(colours = magma(500), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"_filtereddotplot_scaledatabluered.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100), breaks=seq(25,100,5)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()


pdf(file=paste0(output,"_filtereddotplot_scaledatagrayblue.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("gray","blue"), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()
pdf(file=paste0(output,"_filtereddotplot_scaledatainferno.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = inferno(500), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()
head(data)
pdf(file=paste0(output,"_rawplot.pdf"), width=16)
data %>% 
  filter(pct.exp > 5) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp, size=pct.exp)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
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