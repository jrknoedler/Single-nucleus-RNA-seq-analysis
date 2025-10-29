#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Main6_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")

DefaultAssy(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mylevels <- c(6,17,22,34,0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33)

glut <- VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Esr1vlndata.csv"))
pdf(file=paste0(output,"Esr1vln.pdf"), width=30)
ggplot(glutdat, aes(y=Esr1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank())
dev.off()

glut <- VlnPlot(mySeurat, features=c("Gad1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Gad1vlndata.csv"))
pdf(file=paste0(output,"Gad1vln.pdf"), width=30)
ggplot(glutdat, aes(y=Gad1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Cyp19a1vlndata.csv"))
pdf(file=paste0(output,"Cyp19a1vln.pdf"), width=30)
ggplot(glutdat, aes(y=Cyp19a1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Slc17a6"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Slc17a6vlndata.csv"))
pdf(file=paste0(output,"Slc17a6vln.pdf"), width=30)
ggplot(glutdat, aes(y=Slc17a6, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Tac1vlndata.csv"))
pdf(file=paste0(output,"Tac1vln.pdf"), width=30)
ggplot(glutdat, aes(y=Tac1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()