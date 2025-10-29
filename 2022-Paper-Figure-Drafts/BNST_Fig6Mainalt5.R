#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Fig6Revised"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)
library(patchwork)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)

names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)


mainplotgenes <- c("Tac1","Etv1","Maml2","B3gnt2","Nfib","Moxd1","Rreb1","Syne2")
mainplotgenes <- unlist(mainplotgenes)
counts <- mySeurat@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,mainplotgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
BNSTout <- cbind(mainplotgenes,max.df)
BNSTout$CPM <- "CPM"

CPM <- ggplot(BNSTout, aes(y=factor(mainplotgenes, levels=rev(mainplotgenes)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2.2,3.7)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


pdf(file=paste0(output,"Tacr1_peptides.pdf"), height=20, width=25)
VlnPlot(mySeurat, features=c("Tacr1","Cartpt","Cck","Sst","Vip","Tac2"), ncol=1, pt.size=0)
dev.off()

pdf(file=paste0(output,"Sharedmarkers.pdf"), height=60, width=25)
VlnPlot(mySeurat, features=c("Tac1","Etv1","Maml2","B3gnt2","Nfib","Moxd1","Rreb1","Syne2"), ncol=1, pt.size=0)
dev.off()

glut <- VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
glutdat <- glut$data

v1 <- ggplot(glutdat, aes(y=Tac1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))

glut <- VlnPlot(mySeurat, features=c("Etv1"), pt.size=0)
glutdat <- glut$data


v2 <- ggplot(glutdat, aes(y=Etv1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))


glut <- VlnPlot(mySeurat, features=c("Maml2"), pt.size=0)
glutdat <- glut$data


v3 <- ggplot(glutdat, aes(y=Maml2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))



glut <- VlnPlot(mySeurat, features=c("B3gnt2"), pt.size=0)
glutdat <- glut$data


v4 <- ggplot(glutdat, aes(y=B3gnt2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))

glut <- VlnPlot(mySeurat, features=c("Nfib"), pt.size=0)
glutdat <- glut$data

v5 <- ggplot(glutdat, aes(y=Nfib, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))

glut <- VlnPlot(mySeurat, features=c("Moxd1"), pt.size=0)
glutdat <- glut$data

v6 <- ggplot(glutdat, aes(y=Moxd1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))

glut <- VlnPlot(mySeurat, features=c("Rreb1"), pt.size=0)
glutdat <- glut$data

v7 <- ggplot(glutdat, aes(y=Rreb1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))


glut <- VlnPlot(mySeurat, features=c("Syne2"), pt.size=0)
glutdat <- glut$data

v8 <- ggplot(glutdat, aes(y=Syne2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))


pdf(file=paste0(output,"Marker_vln_stacked.pdf"), width=67, height=30)
((v1/v2/v3/v4/v5/v6/v7/v8)|CPM) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Tac1"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vt <- ggplot(glutdat, aes(y=Tac1, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
pdf(file=paste0(output,"Tac1_consplit.pdf"))
vt
dev.off()

glut <- VlnPlot(mySeurat, features=c("Etv1"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs1 <- ggplot(glutdat, aes(y=Etv1, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Maml2"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs2 <- ggplot(glutdat, aes(y=Maml2, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("B3gnt2"), idents=c("18"),pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs3 <- ggplot(glutdat, aes(y=B3gnt2, x=ident, idents=c("18"),fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Nfib"), idents=c("18"),pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs4 <- ggplot(glutdat, aes(y=Nfib, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Moxd1"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs5 <- ggplot(glutdat, aes(y=Moxd1, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Rreb1"),idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs6 <- ggplot(glutdat, aes(y=Rreb1, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none") 

glut <- VlnPlot(mySeurat, features=c("Syne2"),idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs7 <- ggplot(glutdat, aes(y=Syne2, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none") 
blank <- ggplot() +theme_void()
pdf(file=paste0(output,"Markers_Tac1onlycol1.pdf"),height=13)
(vs1/vs2/vs3/vs4) + plot_layout(widths=c("1","0.05")) 
dev.off()
pdf(file=paste0(output,"Markers_Tac1onlycol2.pdf"), height=13)
(vs5/vs6/vs7/blank)  + plot_layout(widths=c("1","0.05")) 
dev.off() 





glut <- VlnPlot(mySeurat, features=c("Maml2"), idents=c("18"), pt.size=0, split.by="Hormone")
glutdat <- glut$data
vs2 <- ggplot(glutdat, aes(y=Maml2, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Arid5b"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs1 <- ggplot(glutdat, aes(y=Arid5b, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Sdk1"), idents=c("18"),pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs3 <- ggplot(glutdat, aes(y=Sdk1, x=ident, idents=c("18"),fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Slc44a5"), idents=c("18"),pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs4 <- ggplot(glutdat, aes(y=Slc44a5, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Sox6"), idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs5 <- ggplot(glutdat, aes(y=Sox6, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none")

glut <- VlnPlot(mySeurat, features=c("Rreb1"),idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs6 <- ggplot(glutdat, aes(y=Rreb1, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none") 

glut <- VlnPlot(mySeurat, features=c("Syne2"),idents=c("18"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
vs7 <- ggplot(glutdat, aes(y=Syne2, x=ident, fill=split)) + geom_violin(trim=TRUE, scale="width") + scale_fill_manual(values=c("royalblue3","deeppink1","chartreuse3"))+theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) + theme(legend.position = "none") 
blank <- ggplot() +theme_void()


pdf(file=paste0(output,"Markers_Tac1onlycol3.pdf"), width=40, height=20)
(vs1/vs2/vs3/vs4)
dev.off()
pdf(file=paste0(output,"Markers_Tac1onlycol4.pdf"), width=40, height=20)
(vs5/vs6/vs7/blank)
dev.off()




cons <- FindConservedMarkers(mySeurat, ident.1="18", grouping.var="Hormone", assay="RNA", test.use="MAST")
write.csv(cons, file=paste0(output, "ConservedTac1Markers.csv"))
pdf(file=paste0(output,"Marker_vln_stackedsplit.pdf"), width=120, height=30)
(vs2/vs3/vs4/vs5/vs6/vs7)
dev.off()

pdf(file=paste0(output,"Marker_vln_stackedsplit2.pdf"), width=120, height=30)
(vs1/vs2/vs3/vs4/vs5/vs6/vs7)
dev.off()