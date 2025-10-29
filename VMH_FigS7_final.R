#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_FigS7_final"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)
library(tidyr)
library(tidyverse)
library(patchwork)
mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
pdf(file=paste0(output,"cckarfeaturesexplit.pdf"), width=12)
FeaturePlot(mySeurat, features=c("Cckar"), cols=c("light gray","red"), split.by="sex")
dev.off()
glut <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev1 <- ggplot(glutdat, aes(y=Cckar, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Mical2"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev2 <- ggplot(glutdat, aes(y=Mical2, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Tmem215"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev3 <- ggplot(glutdat, aes(y=Tmem215, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Mad2l1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev4 <- ggplot(glutdat, aes(y=Mad2l1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Ptp4a1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev5 <- ggplot(glutdat, aes(y=Ptp4a1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Ksr1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data


ev6 <- ggplot(glutdat, aes(y=Ksr1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"FrUPcol1.pdf"), height=17, width=10)
(ev2/ev3/ev4/ev5)
dev.off
pdf(file=paste0(output,"FrUPcol2.pdf"), height=17, width=10)
(ev4/ev5/ev6)
dev.off()
glut <- VlnPlot(mySeurat, features=c("Cgnl1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data
ev7 <- ggplot(glutdat, aes(y=Cgnl1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Popdc3"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev8 <- ggplot(glutdat, aes(y=Popdc3, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Trim36"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev9 <- ggplot(glutdat, aes(y=Trim36, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Klhl1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

ev10 <- ggplot(glutdat, aes(y=Klhl1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
blank <- ggplot() + theme_void()
pdf(file=paste0(output,"FuUPcol1.pdf"), height=17, width=10)
(ev7/ev8/ev9/ev10)
dev.off()
pdf(file=paste0(output,"FuUPcol21.pdf"), height=17, width=10)
(ev9/ev10/blank)
dev.off()





pdf(file=paste0(output,"Sharedmarkers.pdf"), height=25, width=25)
VlnPlot(mySeurat, idents=c("6","20","25"), features=c("Igsf11","Abtb2","Kcnh8","Chsy3","Mpped1","Atp8b1","Zmiz1","Osbpl3","Alk","Trpc3","Slc24a4","Pou6f2"), ncol=4, pt.size=0)
dev.off()

glut <- VlnPlot(mySeurat, features=c("Igsf11","Abtb2","Kcnh8","Chys8","Mpped1","Atp8b1","Zmiz1","Osbpl3","Alk","Trpc3","Slc24a4","Pou6f2"), idents=c("6","20","25"), ncol=1, pt.size=0)


glut <- VlnPlot(mySeurat, features=c("Mpped1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv1 <- ggplot(glutdat, aes(y=Mpped1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Atp8b1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv2 <- ggplot(glutdat, aes(y=Atp8b1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Zmiz1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv3 <- ggplot(glutdat, aes(y=Zmiz1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Nos1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv4 <- ggplot(glutdat, aes(y=Nos1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Tmem163"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv5 <- ggplot(glutdat, aes(y=Tmem163, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Ptprm"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv6 <- ggplot(glutdat, aes(y=Ptprm, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Adam12"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv7 <- ggplot(glutdat, aes(y=Adam12, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Bcl11a"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sfv8 <- ggplot(glutdat, aes(y=Bcl11a, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))



pdf(file=paste0(output,"Sharedfinal_col1.pdf"), height=17, width=10)
(sfv1/sfv2/sfv3/sfv4)
dev.off()

pdf(file=paste0(output,"Sharedfinal_col2.pdf"), height=17, width=10)
(sfv5/sfv6/sfv7/sfv8)
dev.off()







glutdat <- glut$data

sv1 <- ggplot(glutdat, aes(y=Igsf11, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Igsf11"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv1 <- ggplot(glutdat, aes(y=Igsf11, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Abtb2"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv2 <- ggplot(glutdat, aes(y=Abtb2, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Kcnh8"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv3 <- ggplot(glutdat, aes(y=Kcnh8, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Chsy3"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv4 <- ggplot(glutdat, aes(y=Chsy3, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Chsy3"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv5 <- ggplot(glutdat, aes(y=Chsy3, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Mpped1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv6 <- ggplot(glutdat, aes(y=Mpped1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Atp8b1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv7 <- ggplot(glutdat, aes(y=Atp8b1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Zmiz1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv8 <- ggplot(glutdat, aes(y=Zmiz1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Zmiz1"), idents=c("6","20","25"), pt.size=0)

glutdat <- glut$data

sv9 <- ggplot(glutdat, aes(y=Zmiz1, x=ident, fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=30, color="black"), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


pdf(file=paste0(output,"CckarMyo1hmarkers.pdf"), height=200, width=25)
VlnPlot(mySeurat, features=c("Chsy3","Trpc5","Alk","Grk3","Trpc3","Synpr","Ppargc1a","Rreb1","Sema6d","Gfod1","Lin28b","Vav3","Fmn1","Crtac1","Reln","Setdb2","Trpc7","Igsf11","Nrgn","Slc24a4","B3gat2","Ncald","Pou6f2","Chst9","Scml4","Kiz","Stk32b","Cab39l"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"Abtb2Myo1hmarkers.pdf"), height=30, width=25)
VlnPlot(mySeurat, features=c("Slit2","Kcnh8","Abtb2","Megf11","Rxrg","Tox","Osbpl3","Nr2f2"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"unionmarkers.pdf"), height=30, width=25)
VlnPlot(mySeurat, features=c("Ptprm","Atp8b1","Pcsk5","Adam12","Zmiz1","Tmem163","Bcl11a","Mpped1","Gng12","Mcc","Nos1"), pt.size=0)
dev.off()
glut <- VlnPlot(mySeurat, features=c("Mpped1"), pt.size=0)
glutdat <- glut$data


mv1 <- ggplot(glutdat, aes(y=Mpped1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Abtb2"), pt.size=0)
glutdat <- glut$data


mv2 <- ggplot(glutdat, aes(y=Abtb2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,3,1))


glut <- VlnPlot(mySeurat, features=c("Kcnh8"), pt.size=0)
glutdat <- glut$data

mv3 <- ggplot(glutdat, aes(y=Kcnh8, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Trpc3"), pt.size=0)
glutdat <- glut$data

mv4 <- ggplot(glutdat, aes(y=Trpc3, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Synpr"), pt.size=0)
glutdat <- glut$data

mv5 <- ggplot(glutdat, aes(y=Synpr, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Igsf11"), pt.size=0)
glutdat <- glut$data

mv6 <- ggplot(glutdat, aes(Igsf11, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Chsy3"), pt.size=0)
glutdat <- glut$data

mv7 <- ggplot(glutdat, aes(y=Chsy3, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))

#glut <- VlnPlot(mySeurat, features=c("Chsy3"), pt.size=0)
#glutdat <- glut$data

#mv8 <- ggplot(glutdat, aes(y=Rxrg, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank() ,axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"Marker_vln_stacked_female.pdf"), width=80, height=30)
(mv1/mv2/mv3/mv4/mv5/mv6/mv7)
dev.off()
glut <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("19","20","21"), pt.size=0)

glutdat <- glut$data

cckartCT <- ggplot(glutdat, aes(y=Cckar, x=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width")+scale_fill_manual(values=c("blue","deeppink1","green")) + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(legend.position="none") 
pdf(file=paste0(output,"CckartCT_Cckar.pdf"))
cckartCT
dev.off()

v1 <- ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.5,1))


glut <- VlnPlot(mySeurat, features=c("Samd9l"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Samd9lfullvlndata.csv"))

v2 <- ggplot(glutdat, aes(y=Samd9l, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.2,1))


glut <- VlnPlot(mySeurat, features=c("Mylk"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Mylkfullvlndata.csv"))

v3 <- ggplot(glutdat, aes(y=Mylk, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.5,1))


glut <- VlnPlot(mySeurat, features=c("Tmem215"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Tmem215fullvlndata.csv"))

v4 <- ggplot(glutdat, aes(y=Tmem215, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,1.5,1))

glut <- VlnPlot(mySeurat, features=c("St3gal1"), pt.size=0)
glutdat <- glut$data
v5 <- ggplot(glutdat, aes(y=St3gal1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2,1))

glut <- VlnPlot(mySeurat, features=c("Ptp4a1"), pt.size=0)
glutdat <- glut$data
v6 <- ggplot(glutdat, aes(y=Ptp4a1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2,1))

glut <- VlnPlot(mySeurat, features=c("Mical2"), pt.size=0)
glutdat <- glut$data
v7 <- ggplot(glutdat, aes(y=Mical2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,3,1))

glut <- VlnPlot(mySeurat, features=c("Gfra1"), pt.size=0)
glutdat <- glut$data
v8 <- ggplot(glutdat, aes(y=Gfra1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,3,1))


pdf(file=paste0(output,"eDEGvlnsstacked.pdf"), width=60, height=35)
(v1/v3/v4/v5/v6/v7/v8)

dev.off()
glut <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("20"), split.by="Hormone", pt.size=0)

glutdat <- glut$data

cckartCT <- ggplot(glutdat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, scale="width")+scale_fill_manual(values=c("blue","deeppink1","green")) + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=80, color="black"), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))
pdf(file=paste0(output,"CckartCT_Cckar.pdf")
cckartCT
dev.off()

glut <- VlnPlot(mySeurat, features=c("Rprm"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Rprmvlndata.csv"))
pdf(file=paste0(output,"Rprmvln.pdf"), width=40)
ggplot(glutdat, aes(y=Rprm, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()


glut <- VlnPlot(mySeurat, features=c("Iigp1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Iigp1vlndata.csv"))
pdf(file=paste0(output,"Iigp1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Iigp1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Phf21bvlndata.csv"))
pdf(file=paste0(output,"Phf21bvln.pdf"), width=40)
ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Zfp804a"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Zfp804avlndata.csv"))
pdf(file=paste0(output,"Zfp804avln.pdf"), width=40)
ggplot(glutdat, aes(y=Zfp804a, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()



glut <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"CZIcckarvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "cckarCZIvlndata.csv"))
pdf(file=paste0(output,"CckarCZIvln.pdf"))
ggplot(glutdat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()


glut <- VlnPlot(mySeurat, features=c("Myo1h"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"CZImyo1hvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "myo1hCZIvlndata.csv"))
pdf(file=paste0(output,"Myo1hCZIvln.pdf"))
ggplot(glutdat, aes(y=Myo1h, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Phf21b"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"Phf21bCZIvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Phf21bCZIvlndata.csv"))
pdf(file=paste0(output,"Phf21bCZIvln.pdf"))
ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

AvgAllDEGs <- AverageExpression(mySeurat, features=c("Cckar","Rprm","Rnf19a","Iigp1","Maml2","Mc4r","Mad2l1","Tmod1"), assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/Cckarmarkerssavg.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/Cckarmarkerssavg.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
write.csv(data, file=paste0(output,"ScaledCckarMarkers.csv"))
idents <- mySeurat@active.ident
idents <- as.matrix(idents)
head(idents)
pdf(file=paste0(output,"Cckarheatmaptrowclustonly.pdf"))
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Cckarheatmapnoclust.pdf"))
pheatmap(data, cluster_cols=FALSE, cluster_rows=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Cckarheatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()






AvgselectDEGs <- AverageExpression(mySeurat, features=c("Mad2l1","Stat5a","Foxo1","Tmem215","Tnxb","Ezr","Cyp19a1","Gldn","Moxd1","Gal","Pdyn","Il1rap","Rbms3","Cartpt","Calb1","Greb1","Zfp804a","Klhl1","Nxph1","Cdh10","B3galt1","Cab39l","Cgnl1","Popdc3","Trim36","Iqschfp"), assays=c("RNA"))
AvgselectDEGs
write.csv(AvgselectDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/SelectDEGsavg2.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/SelectDEGsavg2.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
pdf(file=paste0(output,"SelectDEGsavgheatmaptrowclustonly.pdf"))
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"SelectDEGsavgheatmapnoclust.pdf"))
pheatmap(data, cluster_cols=FALSE, cluster_rows=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"SelectDEGsavgheatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()