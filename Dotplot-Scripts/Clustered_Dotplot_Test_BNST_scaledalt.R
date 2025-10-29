#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_CelltypesDraft_Malat1regress_Esr1filtered_sexclude_res1.2alt2_peptides"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
#mySeurat <- NormalizeData(mySeurat)
typemarkers <- read.table("DotPlotMarkerLists/BNST_Draft2_malat1regress.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.2)
dev.off()

Pct <- DotPlot(mySeurat, features=c("Cyp19a1","Penk","Pomc","Pdyn","Pnoc","Avp","Oxt","Gast","Cck","Sst","Cst","Rfrp","Npy","Ppy","Prl","Calca","Calcb","Nmb","Grp","Vip","Adcyap1","Ghrh","Crh","Ucn","Tac1","Tac2","Tac3","Nmu","Nms","Nts","Gal","Gnrh","Nps","Nxph1","Nxph2","Nxph4"))
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_filtereddotplot_scaledata.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = viridis(500), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank()) + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"_filtereddotplot_scaledatamagma.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = magma(500), limits=c(-1.5, 2.5))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"_filtereddotplot_scaledatabluered.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
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