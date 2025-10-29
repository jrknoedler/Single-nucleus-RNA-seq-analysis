#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Merge_Fig6storyboard_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_20pcs.rds")
DefaultAssay(mySeurat) <- "RNA"
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
seDEGs <- read.table("Seurat/POA_IndependentAnalysis/POA_Candidate_CelltypSdeg.txt")
unlist(seDEGs)
seDEGs <- as.matrix(seDEGs)
seDEGs.filtered <- intersect(seDEGs, genes.10x)


Pct <- DotPlot(mySeurat, features=seDEGs.filtered)
sink(file="Sinktest2.txt", append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_rawdotplot.pdf"))
data %>% filter(pct.exp > 10) %>% ggplot(aes(x=id, y=features.plot, color=avg.exp, size=pct.exp)) + geom_point() + theme(axis.text.x = element_text(theme(axis.ticks=element_blank()))) + scale_color_gradientn(colors=viridis::viridis)
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