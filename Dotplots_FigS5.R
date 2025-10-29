#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/VMH_IndependentAnalysis/VMH_Merged_Clusterfuck_Marchredo_5pct_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- ScaleData(mySeurat)

tenpct <- read.table("DotPlotMarkerLists/VMH_10pctunique.txt")
tenpct <- unlist(tenpct)
tenpct <- as.matrix(tenpct)
fifteenpct <- read.table("DotPlotMarkerLists/VMH_15pctunique.txt")
fifteenpct <- unlist(fifteenpct)
fifteenpct <- as.matrix(fifteenpct)

pct <- DotPlot(mySeurat, features=c("Gadd45b","Nr5a2","Rhoc","Grm2","Tfcp2l1","Ankrd34c","Chrna3","Fgf2","Ghsr","Hopx","Layn"), idents=c("0","7","10","17","18","19","21","26"))
pctsanity <- DotPlot(mySeurat, features=c("Gadd45b","Nr5a2","Rhoc","Grm2","Tfcp2l1","Ankrd34c","Chrna3","Fgf2","Ghsr","Hopx","Layn"))
data <- pct$data

pdf(file=paste0(output,"_VMH_5pctuctoff_caled.pdf"), width=16)
data %>% 
  filter(pct.exp > 0) %>% ggplot(aes(x=features.plot, y=factor(id), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()
datas <- pctsanity$data
pdf(file=paste0(output,"_VMH_5pctuctoff_scaledsanitycheck.pdf"), width=16)
datas %>% 
  filter(pct.exp > 0) %>% ggplot(aes(x=features.plot, y=factor(id), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pct2 <- DotPlot(mySeurat, features=tenpct, idents=c("2","4","7","9","14","17","19","20","21","23","26","27"))

data2 <- pct2$data
pdf(file=paste0(output,"_VMH_10pctuctoff_caled.pdf"), width=16)
data2 %>% 
  filter(pct.exp > 0) %>% ggplot(aes(x=features.plot, y=factor(id), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

pct3 <-  DotPlot(mySeurat, features=fifteenpct, idents=c("0","4","8","9","14","17","19","20","21","23","24","25","26"))

data3 <- pct3$data
pdf(file=paste0(output,"_VMH_15pctuctoff_caled.pdf"), width=16)
data3 %>% 
  filter(pct.exp > 0) %>% ggplot(aes(x=features.plot, y=factor(id), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()