#!/usr/bin/env Rscript

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)
library(patchwork)

output <- "Seurat/Figure_S5"

BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)
DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)

MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)
DefaultAssay(MeA) <- "RNA"
MeA <- NormalizeData(MeA)

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(VMH)
VMH <- RenameIdents(VMH, new.cluster.ids)
DefaultAssay(VMH) <- "RNA"
VMH <- NormalizeData(VMH)

blank <- ggplot() + theme_void()
v <- VlnPlot(BNST, features=c("Gad1"), idents=c("11"), pt.size=0)
vat <- v$data


VBMG <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

pdf(file=paste0(output,"BNSTgad1.pdf"), height=15)
VBMG
dev.off()
#v <- VlnPlot(POA, features=c("Gad1"), idents=c("13","19","22","24","36"), pt.size=0)
#vat <- v$data
#POAmlevels <- c("13","19","22","24","36")
#VPMG <- ggplot(vat, aes(x=Gad1, y=factor(ident, levels=rev(POAmlevels))))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("13"), pt.size=0)
vat <- v$data
VPMG13 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("22"), pt.size=0)
vat <- v$data
VPMG19 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("22"), pt.size=0)
vat <- v$data
VPMG22 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("24"), pt.size=0)
vat <- v$data
VPMG24 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("36"), pt.size=0)
vat <- v$data
VPMG36 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")




v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("13"), pt.size=0)
vat <- v$data
VPMSL13 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("22"), pt.size=0)
vat <- v$data
VPMSL19 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("22"), pt.size=0)
vat <- v$data
VPMSL22 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("24"), pt.size=0)
vat <- v$data
VPMSL24 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("36"), pt.size=0)
vat <- v$data
VPMSL36 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")




v <- VlnPlot(VMH, features=c("Gad1"), idents=c("17"), pt.size=0)
vat <- v$data

VVMG <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")





v <- VlnPlot(BNST, features=c("Slc17a6"), idents=c("11"), pt.size=0)
vat <- v$data


VBMSL <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("13","19","22","24","36"), pt.size=0)
vat <- v$data
POAmlevels <- c("13","19","22","24","36")
VPMSL <- ggplot(vat, aes(x=Slc17a6, y=factor(ident, levels=rev(POAmlevels))))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Slc17a6"), idents=c("17"), pt.size=0)
vat <- v$data

VVMSL <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")
#gad <- (VBMG/blank/VPMG/blank/VVMG) + plot_layout(heights=c("0.2","0.2","1","0.2","0.2")) & xlim(-0.2,3)
#glut <- (VBMSL/blank/VPMSL/blank/VVMSL) + plot_layout(heights=c("0.2","0.2","1","0.2","0.2")) & xlim(-0.2,3)
#pdf(file=paste0(output,"Male1Vlntest.pdf"))
#(gad|glut)
#dev.off()

v <- VlnPlot(VMH, features=c("Gldn"), idents=c("17"), pt.size=0)
vat <- v$data

VVC17 <- ggplot(vat, aes(x=Gldn, y=factor(ident, levels=rev(POAmlevels))))+ geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(BNST, features=c("Arx"), idents=c("11"), pt.size=0)
vat <- v$data
VBC11 <- ggplot(vat, aes(x=Arx, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Meis2"), idents=c("13"), pt.size=0)
vat <- v$data
VPC13 <- ggplot(vat, aes(x=Meis2, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Thbs4"), idents=c("19"), pt.size=0)
vat <- v$data
VPC19 <- ggplot(vat, aes(x=Thbs4, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Cdh23"), idents=c("22"), pt.size=0)
vat <- v$data
VPC22 <- ggplot(vat, aes(x=Cdh23, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Chst9"), idents=c("24"), pt.size=0)
vat <- v$data
VPC24 <- ggplot(vat, aes(x=Chst9, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Rspo2"), idents=c("36"), pt.size=0)
vat <- v$data
VPC36 <- ggplot(vat, aes(x=Rspo2, y=ident))+ cowplot::theme_nothing() + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


gadsplit <- (VBMG/blank/VPMG13/VPMG19/VPMG22/VPMG24/VPMG36/blank/VVMG) & xlim(-0.2,3)
glutsplit <- (VBMSL/blank/VPMSL13/VPMSL19/VPMSL22/VPMSL24/VPMSL36/blank/VVMSL) & xlim(-0.2,3)
markers <- (VBC11/blank/VPC13/VPC19/VPC22/VPC24/VPC36/blank/VVC17) 
pdf(file=paste0(output,"Male1Vlnfull.pdf"))
(gadsplit|glutsplit|blank|markers) + plot_layout(widths=c("0.8","0.8","1","0.8"))
dev.off()

POAflevels <- c("21","29")
v <- VlnPlot(POA, features=c("Gad1"), idents=c("21"), pt.size=0)
vat <- v$data
VPFG21 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Gad1"), idents=c("29"), pt.size=0)
vat <- v$data
VPFG29 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Gad1"), idents=c("23"), pt.size=0)
vat <- v$data
VMFG <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("21"), pt.size=0)
vat <- v$data
VPFSL21 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Slc17a6"), idents=c("29"), pt.size=0)
vat <- v$data
VPFSL29 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Slc17a6"), idents=c("23"), pt.size=0)
vat <- v$data
VMFSLG <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(MeA, features=c("Antxr2"), idents=c("23"), pt.size=0)
vat <- v$data
VMFM <- ggplot(vat, aes(x=Antxr2, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Kiss1"), idents=c("29"), pt.size=0)
vat <- v$data
VPFM29 <- ggplot(vat, aes(x=Kiss1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(POA, features=c("Npy2r"), idents=c("21"), pt.size=0)
vat <- v$data
VPFM21 <- ggplot(vat, aes(x=Npy2r, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

fgad <- (VMFG/blank/VPFG21/VPFG29) & xlim(-0.2,3)
fglut <- (VMFSLG/blank/VPFSL21/VPFSL29) & xlim(-0.2,3)
fmarkers <- (VMFM/blank/VPFM21/VPFM29)
pdf(file=paste0(output,"Femalevln.pdf"))
(fgad|fglut|blank|fmarkers) + plot_layout(widths=c("0.8","0.8","1","0.8"))
dev.off()


v <- VlnPlot(VMH, features=c("Gad1"), idents=c("6"), pt.size=0)
vat <- v$data
VVFG6 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Gad1"), idents=c("20"), pt.size=0)
vat <- v$data
VVFG20 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Gad1"), idents=c("25"), pt.size=0)
vat <- v$data
VVFG25 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(VMH, features=c("Slc17a6"), idents=c("6"), pt.size=0)
vat <- v$data
VVFSL6 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Slc17a6"), idents=c("20"), pt.size=0)
vat <- v$data
VVFSL20 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Slc17a6"), idents=c("25"), pt.size=0)
vat <- v$data
VVFSL25 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(VMH, features=c("Abtb2"), idents=c("6"), pt.size=0)
vat <- v$data
VVFM6 <- ggplot(vat, aes(x=Abtb2, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Cckar"), idents=c("20"), pt.size=0)
vat <- v$data
VVFM20 <- ggplot(vat, aes(x=Cckar, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(VMH, features=c("Myo1h"), idents=c("25"), pt.size=0)
vat <- v$data
VVFM25 <- ggplot(vat, aes(x=Myo1h, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(MeA, features=c("Gad1"), idents=c("7"), pt.size=0)
vat <- v$data
VMFG7 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Gad1"), idents=c("8"), pt.size=0)
vat <- v$data
VMFG8 <- ggplot(vat, aes(x=Gad1, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Gad1"), idents=c("24"), pt.size=0)
vat <- v$data
VMFG24 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Gad1"), idents=c("32"), pt.size=0)
vat <- v$data
VMFG32 <- ggplot(vat, aes(x=Gad1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(MeA, features=c("Slc17a6"), idents=c("7"), pt.size=0)
vat <- v$data
VMFSL7 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Slc17a6"), idents=c("8"), pt.size=0)
vat <- v$data
VMFSL8 <- ggplot(vat, aes(x=Slc17a6, y=ident))+ geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Slc17a6"), idents=c("24"), pt.size=0)
vat <- v$data
VMFSL24 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Slc17a6"), idents=c("32"), pt.size=0)
vat <- v$data
VMFSL32 <- ggplot(vat, aes(x=Slc17a6, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="dark gray") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(MeA, features=c("Npr3"), idents=c("7"), pt.size=0)
vat <- v$data
VMFM7 <- ggplot(vat, aes(x=Npr3, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Isl1"), idents=c("8"), pt.size=0)
vat <- v$data
VMFM8 <- ggplot(vat, aes(x=Isl1, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Ebf3"), idents=c("24"), pt.size=0)
vat <- v$data
VMFM24 <- ggplot(vat, aes(x=Ebf3, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(MeA, features=c("Adgrg2"), idents=c("32"), pt.size=0)
vat <- v$data
VMFM32 <- ggplot(vat, aes(x=Adgrg2, y=ident)) + geom_violin(trim=TRUE, scale="width", fill="black") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

emarkers <- (VMFM7/VMFM8/VMFM24/VMFM32/blank/VVFM6/VVFM20/VVFM25)
eglut <- (VMFSL7/VMFSL8/VMFSL24/VMFSL32/blank/VVFSL6/VVFSL20/VVFSL25) & xlim(-0.2,3)
egad <- (VMFG7/VMFG8/VMFG24/VMFG32/blank/VVFG6/VVFG20/VVFG25) & xlim(-0.2,3)
pdf(file=paste0(output,"estrusvln.pdf"))
(egad|eglut|blank|emarkers)  + plot_layout(widths=c("0.8","0.8","1","0.8"))
dev.off()
