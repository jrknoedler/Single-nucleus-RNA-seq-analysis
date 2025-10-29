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


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test.rds")
DefaultAssay(mySeurat) <- "RNA"
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
seDEGs <- read.table("Seurat/POA_IndependentAnalysis/POA_Candidate_CelltypSdeg.txt")
unlist(seDEGs)
seDEGs <- as.matrix(seDEGs)
seDEGs.filtered <- intersect(seDEGs, genes.10x)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Pct <- DotPlot(mySeurat, features=seDEGs.filtered)
sink(file="Sinktest2.txt", append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)

pdf(file=paste0(output,"_rawdotplot.pdf"), width=12)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=id, color=avg.exp, size=pct.exp)) + 
  geom_point() + 
  scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

# make data square to calculate euclidean distance
mat <- data %>% select(-avg.exp.scaled, -pct.exp) %>%
	pivot_wider(names_from=id, values_from=avg.exp) %>%
	data.frame() # make df as tibbles -> matrix annoying`
row.names(mat) <- mat$features.plot
mat <- mat[,-1]
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix

pdf(file=paste0(output,"clusteredDotplot.pdf"))
dotplot <- data  %>% 
  mutate(id = factor(id, levels = clust$order)) %>% 
  filter(pct.exp > 25) %>% 
  ggplot(aes(x=features.plot, y=id, color=avg.exp, size=pct.exp)) + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
dotplot
dev.off()