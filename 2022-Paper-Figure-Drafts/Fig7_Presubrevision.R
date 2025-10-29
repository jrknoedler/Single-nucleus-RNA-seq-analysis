#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "X:/Seurat/VMH_IndependentAnalysis/VMH_FigS7Main7_final_revised"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)
library(tidyr)
library(tidyverse)
library(patchwork)
library(ggnewscale)
mySeurat <- readRDS("X:/Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")

DefaultAssay(mySeurat) <- "RNA"





new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
mySeurat <- NormalizeData(mySeurat)

Male <- subset(mySeurat, subset=Hormone=="Intact")
Primed <- subset(mySeurat, subset=Hormone=="Primed")
Unprimed <- subset(mySeurat, subset=Hormone=="Unprimed")



v <- VlnPlot(Male, features=c("Cckar"), pt.size=0)
vat <- v$data


vve1 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("dark gray",1), rep("dark gray",7)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=0.5), axis.ticks=element_blank(), axis.line.y=element_blank()) + scale_y_continuous(breaks=seq(0,4,1))+theme(legend.position="none") + scale_x_discrete(drop=FALSE)

v <- VlnPlot(Primed, features=c("Cckar"), pt.size=0)
vat <- v$data


vve2 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",7)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=0.5), axis.ticks=element_blank(), axis.line.y=element_blank()) + scale_y_continuous(breaks=seq(0,4,1))+theme(legend.position="none") + scale_x_discrete(drop=FALSE)

v <- VlnPlot(Unprimed, features=c("Cckar"), pt.size=0)
vat <- v$data


vve3 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("dark gray",1), rep("dark gray",7)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=0.5), axis.ticks=element_blank(), axis.line.y=element_blank()) + scale_y_continuous(breaks=seq(0,4,1))+theme(legend.position="none") + scale_x_discrete(drop=FALSE)

pdf(file=paste0(output,"Cckarcondplitstaked.pdf"), width=20, height=5)
(vve1/vve2/vve3) & theme(axis.ticks=element_blank())
dev.off()
