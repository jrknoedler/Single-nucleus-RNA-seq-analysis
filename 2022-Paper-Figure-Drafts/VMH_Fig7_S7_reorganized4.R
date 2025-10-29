#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_FigS7Main7_final"

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
library(ggnewscale)
mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
blank <- ggplot() + theme_void()
DefaultAssay(mySeurat) <- "RNA"

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
mySeurat <- NormalizeData(mySeurat)


markerscpm <- c("Cckar","Trim36")
markerscpm <- unlist(markerscpm)

edegscpm <- c("Arg2","Tmem215","Gfra1","Mylk","Nfil3","Tmod1","Phf21b")
edegscpm <- unlist(edegscpm)

commoncpm <- c("Mpped1","Atp8b1","Trpc3","Synpr","Igsf11","Chsy3")
commoncpm <- unlist(commoncpm)

counts <- mySeurat@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)
head(CPM)
genes.m <- CPM[,markerscpm]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
markerscpmout <- cbind(markerscpm,max.df)
markerscpmout$CPM <- "CPM"

genes.m <- CPM[,edegscpm]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
edegscpmout <- cbind(edegscpm,max.df)
edegscpmout$CPM <- "CPM"

genes.m <- CPM[,commoncpm]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
commoncpmout <- cbind(commoncpm ,max.df)
commoncpmout$CPM <- "CPM"

#pdf(file=paste0(output,"DEGcpmscalebar.pdf"), height=40)
#ggplot(Mainout, aes(y=factor(mainplotgenes, levels=rev(mainplotgenes)), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
#dev.off()

markerCPMplot <- ggplot(markerscpmout, aes(x=factor(markerscpm, levels=markerscpm), y=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2,3.3)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank())  + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
edegCPMplot <- ggplot(edegscpmout, aes(x=factor(edegscpm, levels = edegscpm), y=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2,3.3)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank()) 
commonCPMplot <- ggplot(commoncpmout, aes(x=factor(commoncpm, levels=commoncpm), y=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2.5,3.7)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file=paste0(output,"markerCMPscalebar.pdf"))
ggplot(markerscpmout, aes(x=factor(markerscpm, levels=markerscpm), y=CPM, fill=log10(max))) +cowplot::theme_cowplot() + geom_point(size=30, shape=21) + scale_fill_gradientn(colors=c("white","red")) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank())  + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"edegsCPMscalebar.pdf"))
ggplot(edegscpmout, aes(x=factor(edegscpm, levels = edegscpm), y=CPM, fill=log10(max))) +cowplot::theme_cowplot() + geom_point(size=30, shape=21) + scale_fill_gradientn(colors=c("white","red")) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank()) 
dev.off()

pdf(file=paste0(output,"commonCPMscalebar.pdf"))
ggplot(commoncpmout, aes(x=factor(commoncpm, levels=commoncpm), y=CPM, fill=log10(max))) +cowplot::theme_cowplot() + geom_point(size=30, shape=21) + scale_fill_gradientn(colors=c("white","red")) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

pdf(file=paste0(output,"doubleplotMyo1h.pdf"),width=18)

FeaturePlot(mySeurat, features=c("Cckar","Myo1h"), blend=TRUE, order=TRUE, cols=c("light grey","dark red","dark green"))
dev.off()

pdf(file=paste0(output,"doubleplotTrim36.pdf"),width=18)

FeaturePlot(mySeurat, features=c("Cckar","Trim36"), blend=TRUE, order=TRUE, cols=c("light grey","dark red","dark green"))
dev.off()


pdf(file=paste0(output,"doubleplotMyo1hsexsplit.pdf"),width=18, height=12)

FeaturePlot(mySeurat, features=c("Cckar","Myo1h"), blend=TRUE, order=TRUE, split.by="sex",cols=c("light grey","dark red","dark green"))
dev.off()

pdf(file=paste0(output,"doubleplotTrim36sexsplit.pdf"),width=18, height=12)

FeaturePlot(mySeurat, features=c("Cckar","Trim36"), blend=TRUE, order=TRUE, split.by="sex",cols=c("light grey","dark red","dark green"))
dev.off()

t1 <- FeaturePlot(mySeurat, features=c("Cckar","Trim36"), blend=TRUE, order=TRUE, split.by="sex",cols=c("light grey","dark red","dark green"))
head(t1)
#saveRDS(t1, file="trim36fullsplit.rds")

t2 <- FeaturePlot(mySeurat, features=c("Cckar","Myo1h"), blend=TRUE, order=TRUE, split.by="sex",cols=c("light grey","dark red","dark green"))
#saveRDS(t2, file="myo1hfullsplit.rds")

edegmarkers <- read.table("Genelists/Cckar_eDEGmarkers.txt", header=FALSE)
edegmarkers <- unlist(edegmarkers)
pdf(file=paste0(output,"Tonsofvagplots.pdf"), width=20,height=300)
VlnPlot(mySeurat, features=edegmarkers,ncol=1,pt.size=0)
dev.off()

commonmarkers <- read.table("Genelists/CckarMyo1hcommon.txt", header=FALSE)
common <- unlist(commonmarkers)
pdf(file=paste0(output,"evenmoreofvagplots.pdf"), width=20,height=300)
VlnPlot(mySeurat, features=common,ncol=1,pt.size=0)
dev.off()

Male <- subset(mySeurat, subset=sex=="Male")
Female <- subset(mySeurat, subset=sex=="Female")
Primed <- subset(mySeurat, subset=Hormone=="Primed")
Unprimed <- subset(mySeurat, subset=Hormone=="Unprimed")

Pc <- FeaturePlot(Primed, features=c("Cckar"))
Uc <- FeaturePlot(Unprimed, features=c("Trim36"))
Pcdat <- Pc$data
Ucdat <- Uc$data

Pcdat_1 <- Pcdat[Pcdat$Cckar == 0,]
Pcdat_2 <- Pcdat[Pcdat$Cckar >0,]

Ucdat_1 <- Ucdat[Ucdat$Trim36 == 0,]
Ucdat_2 <- Ucdat[Ucdat$Trim36 >0,]

cuf <- ggplot() + geom_point(data=Pcdat_1, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) + geom_point(data=Pcdat_2, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) +scale_colour_gradientn(colours = c("gray","deeppink1"), limits=c(0, 1.2))+ cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())  
tot <- cuf + new_scale_color() +geom_point(data=Ucdat_1, aes(x=UMAP_1, y=UMAP_2, color=Trim36), size=0.3) + geom_point(data=Ucdat_2, aes(x=UMAP_1, y=UMAP_2, color=Trim36),size=0.3) + scale_color_gradientn(colors=c("gray","green"), limits=c(0, 1.6)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())

Mcc <- FeaturePlot(Male, features=c("Cckar"))
Mmc <- FeaturePlot(Male, features=c("Trim36"))
Mccdat <- Mcc$data
Mmcdat <- Mmc$data

Mccdat_1 <- Mccdat[Mccdat$Cckar == 0,]
Mccdat_2 <- Mccdat[Mccdat$Cckar >0,]

Mmcdat_1 <- Mmcdat[Mmcdat$Trim36 == 0,]
Mmcdat_2 <- Mmcdat[Mmcdat$Trim36 >0,]
Mcuf <- ggplot() + geom_point(data=Mccdat_1, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) + geom_point(data=Mccdat_2, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) +scale_colour_gradientn(colours = c("gray","deeppink1"), limits=c(0, 1.2))+ cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())  +theme(legend.position="none")
Mtot <- Mcuf + new_scale_color() +geom_point(data=Mmcdat_1, aes(x=UMAP_1, y=UMAP_2, color=Trim36), alpha=0.5, size=0.3) + geom_point(data=Mmcdat_2, aes(x=UMAP_1, y=UMAP_2, color=Trim36), alpha=0.5,size=0.3) + scale_color_gradientn(colors=c("gray","green"), limits=c(0, 1.6)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())+theme(legend.position="none")



pdf(file=paste0(output,"fstacktest.pdf"))
tot
dev.off()

pdf(file=paste0(output,"fullUMAPtest.pdf"), width=15)
(Mtot|blank|tot) + plot_layout(widths=c(1,0.1,1))
dev.off()

Mc <- FeaturePlot(Male, features=c("Cckar"), order=TRUE)
mcdat <- Mc$data
head(mcdat)
mcdat_1 <- mcdat[mcdat$Cckar == 0,]
mcdat_2 <- mcdat[mcdat$Cckar >0,]
Fc <- FeaturePlot(Female, features=c("Cckar"))
fcdat <- Fc$data
fcdat_1 <- fcdat[fcdat$Cckar == 0,]
fcdat_2 <- fcdat[fcdat$Cckar > 0,]

dim <- DimPlot(mySeurat, reduction="umap")
dimdat <- dim$data
idu <- ggplot() + geom_point(data=dimdat, aes(x=UMAP_1, y=UMAP_2, color=ident), size=0.3) + scale_color_manual(values=c(rep("dark gray", 5), rep("deeppink1",1), rep("dark gray",13), rep("orange",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))+cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank()) +theme(legend.position="none") 

head(dimdat)
po1 <- ggplot() + geom_point(data=mcdat_1, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) + geom_point(data=mcdat_2, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) +scale_colour_gradientn(colours = c("light blue","dark red"), limits=c(0, 1.2))+ cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())  
po2 <- ggplot() + geom_point(data=fcdat_1, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) + geom_point(data=fcdat_2, aes(x=UMAP_1, y=UMAP_2, color=Cckar), size=0.3) +scale_colour_gradientn(colours = c("light blue","dark red"), limits=c(0, 1.2))+ cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank())+theme(axis.title=element_blank(), axis.ticks=element_blank()) +theme(legend.position="none")

col1 <- (blank/idu/blank) + plot_layout(heights=c(0.275,0.45,0.275))
col2 <- (blank/blank/blank) 
col3 <- po1/blank/po2 + plot_layout(heights=c(0.45,0.1,0.45))
pdf(file=paste0(output,"UMAP_custordertest.pdf"), width=8)
col1|col3 
dev.off()
p1 <- ggplot(mcdat, aes(x=UMAP_1, y=UMAP_2, color=Cckar)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light gray","red"), limits=c(0, 1)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank())


p2 <- ggplot(fcdat, aes(x=UMAP_1, y=UMAP_2, color=Cckar)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light gray","red"), limits=c(0, 1)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank())

pdf(file=paste0(output,"UmapCckar.pdf"), width=12)
p1 + p2
dev.off()



v <- VlnPlot(mySeurat, features=c("Cckar"), pt.size=0)
vat <- v$data


vv1 <- ggplot(vat, aes(x=Cckar, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))+ scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

pdf(file=paste0(output,"vertvlntest.pdf"), height=15)
vv1
dev.off()
vlnlevels <- c(6,20,25)

v <- VlnPlot(mySeurat, features=c("Abtb2"), pt.size=0)
vat <- v$data



vv2 <- ggplot(vat, aes(x=Abtb2, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Trim36"), pt.size=0)
vat <- v$data


vv3 <- ggplot(vat, aes(x=Trim36, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

pdf(file=paste0(output,"vertvln.pdf"), height=15)
vv1
dev.off()


pdf(file=paste0(output,"vertvlntest.pdf"), height=30, width=9.7)
((vv1|vv3)/ markerCPMplot) + plot_layout(heights=c(1,0.05))
dev.off()

markers <- (vv1|vv3)

v <- VlnPlot(mySeurat, features=c("Arg2"), pt.size=0)
vat <- v$data


vve1 <- ggplot(vat, aes(x=Arg2, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Tmem215"), pt.size=0)
vat <- v$data


vve2 <- ggplot(vat, aes(x=Tmem215, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Gfra1"), pt.size=0)
vat <- v$data


vve3 <- ggplot(vat, aes(x=Gfra1, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")



v <- VlnPlot(mySeurat, features=c("Mylk"), pt.size=0)
vat <- v$data


vve4 <- ggplot(vat, aes(x=Mylk, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Nfil3"), pt.size=0)
vat <- v$data


vve5 <- ggplot(vat, aes(x=Nfil3, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Tmod1"), pt.size=0)
vat <- v$data


vve6 <- ggplot(vat, aes(x=Tmod1, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Tmod1"), pt.size=0)
vat <- v$data


vve7 <- ggplot(vat, aes(x=Phf21b, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
vat <- v$data


vve7 <- ggplot(vat, aes(x=Phf21b, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


pdf(file=paste0(output,"M7_eDEGsrow.pdf"), width=34, height=30)
((vve1|vve2|vve3|vve4|vve5|vve6|vve7) /edegCPMplot) + plot_layout(heights=c(1,0.05))
dev.off()

eDEGs <- (vve1|vve2|vve3|vve4|vve5|vve6|vve7)

pdf(file=paste0(output,"M7AB_Combbined.pdf"), width=16, height=30)
markers + eDEGs 
dev.off()

v <- VlnPlot(mySeurat, features=c("Mpped1"), pt.size=0)
vat <- v$data


vvm1 <- ggplot(vat, aes(x=Mpped1, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Atp8b1"), pt.size=0)
vat <- v$data


vvm2 <- ggplot(vat, aes(x=Atp8b1, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Kcnh8"), pt.size=0)
vat <- v$data


vvm3 <- ggplot(vat, aes(x=Kcnh8, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Trpc3"), pt.size=0)
vat <- v$data


vvm4 <- ggplot(vat, aes(x=Trpc3, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Synpr"), pt.size=0)
vat <- v$data


vvm5 <- ggplot(vat, aes(x=Synpr, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Igsf11"), pt.size=0)
vat <- v$data


vvm6 <- ggplot(vat, aes(x=Igsf11, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Chsy3"), pt.size=0)
vat <- v$data


vvm7 <- ggplot(vat, aes(x=Chsy3, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Bcl11a"), pt.size=0)
vat <- v$data


vvm8 <- ggplot(vat, aes(x=Bcl11a, y=factor(ident, levels=rev(mylevels)), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("dark gray",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")


v <- VlnPlot(mySeurat, features=c("Cckar"), pt.size=0)
vat <- v$data


pdf(file=paste0(output,"S7_fmarkersrow.pdf"), width=40, height=30)
((vvm1|vvm2|vvm4|vvm5|vvm6|vvm7)/commonCPMplot) + plot_layout(heights=c(1,0.05))
dev.off()

fv1 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Arg2"), pt.size=0)
vat <- v$data


fv2 <- ggplot(vat, aes(y=Arg2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmem215"), pt.size=0)
vat <- v$data


fv3 <- ggplot(vat, aes(y=Tmem215, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Gfra1"), pt.size=0)
vat <- v$data


fv4 <- ggplot(vat, aes(y=Gfra1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Mylk"), pt.size=0)
vat <- v$data


fv5 <- ggplot(vat, aes(y=Mylk, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
vat <- v$data


fv6 <- ggplot(vat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
v <- VlnPlot(mySeurat, features=c("Nfil3"), pt.size=0)
vat <- v$data


fv7 <- ggplot(vat, aes(y=Nfil3, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmod1"), pt.size=0)
vat <- v$data


fv8 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"S1_eDEGsstacked.pdf"), width=30, height=20)
fv1/fv2/fv3/fv4/fv5/fv6/fv7/fv8
dev.off()


v <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("6","20","25"))

vat <- v$data

vm1 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Abtb2"), idents=c("6","20","25"))

vat <- v$data

vm2 <- ggplot(vat, aes(y=Abtb2, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("6","20","25"))

vat <- v$data

vm3 <- ggplot(vat, aes(y=Cyp19a1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Myo1h"), idents=c("6","20","25"))

vat <- v$data

vm4 <- ggplot(vat, aes(y=Myo1h, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
blank <- ggplot() + theme_void()
pdf(file=paste0(output,"MarkerVlnStack.pdf"), height=13)
vm1/vm2/vm4/blank
dev.off()


v <- VlnPlot(mySeurat, features=c("Arg2"), idents=c("6","20","25"))

vat <- v$data

ve1 <- ggplot(vat, aes(y=Arg2, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Tmem215"), idents=c("6","20","25"))

vat <- v$data

ve2 <- ggplot(vat, aes(y=Tmem215, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Gfra1"), idents=c("6","20","25"))

vat <- v$data

ve3 <- ggplot(vat, aes(y=Gfra1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Mylk"), idents=c("6","20","25"))

vat <- v$data

ve4 <- ggplot(vat, aes(y=Mylk, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
pdf(file=paste0(output,"FrUP1.pdf"), height=13)
ve1/ve2/ve3/ve4
dev.off()

v <- VlnPlot(mySeurat, features=c("Phf21b"), idents=c("6","20","25"))

vat <- v$data

ve5 <- ggplot(vat, aes(y=Phf21b, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Nfil3"), idents=c("6","20","25"))

vat <- v$data

ve6 <- ggplot(vat, aes(y=Nfil3, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
v <- VlnPlot(mySeurat, features=c("Tmod1"), idents=c("6","20","25"))

vat <- v$data

ve7 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmod1"), idents=c("6","20","25"))

vat <- v$data

ve8 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"FrUP2.pdf"), height=13)
ve5/ve6/ve7/blank
dev.off()

v <- VlnPlot(mySeurat, features=c("Cgnl1"), idents=c("6","20","25"))

vat <- v$data

vd1 <- ggplot(vat, aes(y=Cgnl1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Popdc3"), idents=c("6","20","25"))

vat <- v$data

vd2 <- ggplot(vat, aes(y=Popdc3, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Trim36"), idents=c("6","20","25"))

vat <- v$data

vd3 <- ggplot(vat, aes(y=Trim36, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Klhl1"), idents=c("6","20","25"))

vat <- v$data

vd4 <- ggplot(vat, aes(y=Klhl1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"FuUP1.pdf"), width=13)
vd1/vd2/vd3/vd4
dev.off()


counts <- mySeurat@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)
cckar <- CPM[,"Cckar"]
head(cckar)
max(cckar)

genes <- c("Mpped1","Atp8b1","Trpc3","Synpr","Igsf11","Chsy3","Cckar","Myo1h","Arg2","Tmem215","Gfra1","Mylk","Nfil3","Tmod1","Phf21b")
genes
genes.m <- CPM[,genes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
genes <- data.frame(genes)
max.df <- data.frame(max)
out <- cbind(genes,max.df)
out
write.csv(out, file=paste0(output,"CPMcompiled.csv"))
