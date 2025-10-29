#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_PrimedupClusterhighlightmvp_regressed"

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
Maleup <- read.table("Genelists/VMH_PrimedupConsensusFig6.txt", header=FALSE)
Maleup <- unlist(Maleup)
Maleup <- unique(Maleup)
Maleup <- as.matrix(Maleup)
mySeurat <- ScaleData(mySeurat)
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=Maleup, dot.min=0.20)
dev.off()

Pct <- DotPlot(mySeurat, features=Maleup)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
byCckar <- data %>% filter(id==21)   
byCckar
byMaml2 <- data %>%filter(features.plot=="Maml2")
sortedfeature <- byMaml2[order(-byMaml2$avg.exp.scaled),]
sorted <- byCckar[order(-byCckar$pct.exp),]
sorted
pdf(file=paste0(output,"_rawdotplotscaled.pdf"), width=16)
data %>% 
  filter(pct.exp > 20) %>% ggplot(aes(x=(factor(features.plot,levels=sorted$features.plot)), y=factor(id, levels=rev(sortedfeature$id)), color=avg.exp.scaled, size=pct.exp)) + 
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