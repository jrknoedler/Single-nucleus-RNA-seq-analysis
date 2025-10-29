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
library(patchwork)

output <- "Seurat/UniqueDEGDots"

BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
BNST <- ScaleData(BNST)

BNSTlevels <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)
BNSTlevels <- rev(BNSTlevels)

BNSTDEGs <- c("Rbm20","Plekhg1","Prok2","Cck","Phf21b","Slco2a1","Tac1","Asic4","Kank1","Sst","Abca1","Amot","Fam107b","Prkcd","Sertm1","Tns1","Sytl4","Npr3","Plscr4")
Pct <- DotPlot(BNST, features=BNSTDEGs)

data <- Pct$data
BNSTp1 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=BNSTlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(25,100), breaks=seq(25,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') + 
  theme(axis.ticks = element_blank())

BNSTp2 <- data %>% 
  filter(pct.exp > 20) %>% ggplot(aes(x=features.plot, y=factor(id, levels=BNSTlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(20,100), breaks=seq(20,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

BNSTp3 <- data %>% 
  filter(pct.exp > 5) %>% ggplot(aes(x=features.plot, y=factor(id, levels=BNSTlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(5,100), breaks=seq(5,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

BNSTp4 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=BNSTlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(0,100), breaks=seq(0,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

MeA <- ScaleData(MeA)

MeADEGs <- c("Sst","Map3k15","Tspan18","Vgll3","Efemp1")
Pct <- DotPlot(MeA, features=MeADEGs)
MeAlevels <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
MeAlevels <- rev(MeAlevels)

data <- Pct$data
MeAp1 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(25,100), breaks=seq(25,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

MeAp2 <- data %>% 
  filter(pct.exp > 20) %>% ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(20,100), breaks=seq(20,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

MeAp3 <- data %>% 
  filter(pct.exp > 5) %>% ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(5,100), breaks=seq(5,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

MeAp4 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(0,100), breaks=seq(0,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
POA <- ScaleData(POA)

POADEGs <- c("Efemp1","Mob3b","Moxd1","St3gal1","Thbs4","Mrc2","Ttn","Neb","Kiss1","Kl","Slc17a8","Abtb2","Penk")

POAlevels <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39)
POAlevels <- rev(POAlevels)

Pct <- DotPlot(POA, features=POADEGs)

data <- Pct$data


POAp1 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(25,100), breaks=seq(25,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

POAp2 <- data %>% 
  filter(pct.exp > 20) %>% ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(20,100), breaks=seq(20,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

POAp3 <- data %>% 
  filter(pct.exp > 5) %>% ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(5,100), breaks=seq(5,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

POAp4 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(0,100), breaks=seq(0,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
VMH <- ScaleData(VMH)
specDEGs <- c("Cckar","Tnfaip8l3","Mme","Cgnl1","Popdc3","Trim36")
Pct <- DotPlot(VMH, features=specDEGs)
data <- Pct$data

VMHlevels <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)
VMHlevels <- rev(VMHlevels)


VMHp1 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(25,100), breaks=seq(25,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

VMHp2 <- data %>% 
  filter(pct.exp > 20) %>% ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(20,100), breaks=seq(20,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

VMHp3 <- data %>% 
  filter(pct.exp > 5) %>% ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(5,100), breaks=seq(5,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())

VMHp4 <- data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 5))+ scale_size(limits=c(0,100), breaks=seq(0,100,25)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + theme(legend.position="none") +
  ylab('') +
  theme(axis.ticks = element_blank())
blank <- ggplot() + theme_void()
pdf(file=paste0(output,"25filt.pdf"), width=22, height=8.5)
(BNSTp1|blank|MeAp1|blank|POAp1|blank|VMHp1) + plot_layout(widths=c(1,0.1,0.26,0.1,0.68,0.1,0.31))
dev.off()
pdf(file=paste0(output,"20filt.pdf"), width=60, height=20)
(BNSTp2|blank|MeAp2|blank|POAp2|blank|VMHp2) + plot_layout(widths=c(1,0.1,0.26,0.1,0.68,0.1,0.31))
dev.off()
pdf(file=paste0(output,"5filt.pdf"), width=60, height=20)
(BNSTp3|blank|MeAp3|blank|POAp3|blank|VMHp3) + plot_layout(widths=c(1,0.1,0.26,0.1,0.68,0.1,0.31))
dev.off()
pdf(file=paste0(output,"nofilt.pdf"), width=60, height=20)
(BNSTp4|blank|MeAp4|blank|POAp4|blank|VMHp4) + plot_layout(widths=c(1,0.1,0.26,0.1,0.68,0.1,0.31))
dev.off()