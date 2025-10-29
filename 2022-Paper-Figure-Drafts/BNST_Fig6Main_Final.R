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
output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs6dotplot"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

#new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
#mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)

#names(new.cluster.ids) <- levels(mySeurat)
#mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
#head(mySeurat[[]])

mainplotgenes <- c("Tac1","Adra1a","Greb1","Sox5","Synpr","Pura","Socs2","Mical2","Cck","Phf21b","Fam19a2","Kcnh7","Tiparp","Col25a1")
mainplotgenes <- unlist(mainplotgenes)
counts <- mySeurat@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)
head(CPM)
genes.m <- CPM[,mainplotgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
Mainout <- cbind(mainplotgenes,max.df)
Mainout$CPM <- "CPM"

pdf(file=paste0(output,"DEGcpmscalebar.pdf"), height=40)
ggplot(Mainout, aes(y=factor(mainplotgenes, levels=rev(mainplotgenes)), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=30) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25), breaks=c(2,5)+labels=c(2,5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()



col1 <- c("Tac1","Adra1a","Greb1","Sox5")
col1 <- unlist(col1)
col2 <- c("Synpr","Pura","Socs2","Mical2")
col2 <- unlist(col2)
col3 <- c("Cck","Phf21b","Fam19a2","Kcnh7")
col3 <- unlist(col3)
col4 <- c("Tiparp","Col25a1","Esr1","Ar")
col4 <- unlist(col4)

genes.m <- CPM[,col1]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
Col1out <- cbind(col1,max.df)
Col1out$CPM <- "CPM"

genes.m <- CPM[,col2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
Col2out <- cbind(col2,max.df)
Col2out$CPM <- "CPM"

genes.m <- CPM[,col3]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
Col3out <- cbind(col3,max.df)
Col3out$CPM <- "CPM"

genes.m <- CPM[,col4]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
Col4out <- cbind(col4,max.df)
Col4out$CPM <- "CPM"

head(Col1out)
col1CPM <- ggplot(Col1out, aes(y=factor(col1, levels=rev(col1)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=10) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
col2CPM <- ggplot(Col2out, aes(y=factor(col2, levels=rev(col2)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=10) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank()) 
col3CPM <- ggplot(Col3out, aes(y=factor(col3, levels=rev(col3)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=10) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
col4CPM <- ggplot(Col4out, aes(y=factor(col4, levels=rev(col4)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=10) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file=paste0(output,"redscale.pdf"))
col4CPM <- ggplot(Col4out, aes(y=factor(col4, levels=rev(col4)), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=10) + scale_color_gradientn(colors=c("white","red"), limits=c(2.25,4.25), breaks=c(2,5), labels=c(2,5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)


subset <- subset(x=mySeurat, subset=Hormone=="Primed", invert=TRUE)
subset2 <- subset(x=mySeurat, subset=Hormone=="Intact", invert=TRUE)
subset3 <- subset(x=mySeurat, subset=Hormone=="Unprimed", invert=TRUE)

glut <- VlnPlot(mySeurat, features=c("Tac1"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v1 <- ggplot(glutdat, aes(y=Tac1, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
pdf(file=paste0(output,"tac1test.pdf"))
v1
dev.off()

glut <- VlnPlot(mySeurat, features=c("Synpr"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v2 <- ggplot(glutdat, aes(y=Synpr, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Cck"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v3 <- ggplot(glutdat, aes(y=Cck, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Socs2"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v4 <- ggplot(glutdat, aes(y=Socs2, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) 

glut <- VlnPlot(mySeurat, features=c("Adra1a"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v5 <- ggplot(glutdat, aes(y=Adra1a, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Pura"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v6 <- ggplot(glutdat, aes(y=Pura, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1)) 




glut <- VlnPlot(mySeurat, features=c("Greb1"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
v7 <- ggplot(glutdat, aes(y=Greb1, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Mical2"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v8 <- ggplot(glutdat, aes(y=Mical2, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Phf21b"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v9 <- ggplot(glutdat, aes(y=Phf21b, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Kcnh7"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v10 <- ggplot(glutdat, aes(y=Kcnh7, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Sox5"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v11 <- ggplot(glutdat, aes(y=Sox5, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Fam19a2"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data
head(glutdat)

v12 <- ggplot(glutdat, aes(y=Fam19a2, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Etl4"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v13 <- ggplot(glutdat, aes(y=Etl4, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Cnksr1"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v14 <- ggplot(glutdat, aes(y=Cnksr1, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Ap1s2"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v15 <- ggplot(glutdat, aes(y=Ap1s2, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Ap1s2"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v16 <- ggplot(glutdat, aes(y=Ap1s2, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Tiparp"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v17 <- ggplot(glutdat, aes(y=Tiparp, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Col25a1"), idents=c("17"), pt.size=0, split.by="Hormone")

glutdat <- glut$data

v18 <- ggplot(glutdat, aes(y=Col25a1, x=ident, fill=split)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("blue","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

blank <- ggplot() + theme_void()
pdf(file=paste0(output,"Tac1DEGs_col1.pdf"), height=17)
v1/v5/v7/v11/v2/v6/v14 
dev.off()


pdf(file=paste0(output,"Tac1DEGs_col2.pdf"), height=13)
v4/v8/v3/v9/v12/v10 
dev.off()

#pdf(file=paste0(output,"Tac1DEAs_col1b.pdf"),height=17)
#(v2/blank/v6/blank/v4/blank/v8) + plot_layout(heights=c(1,0.4,1,0.4,1,0.4,1))
#dev.off()
col1b <- (blank/v2/blank/v6/blank/v4/blank/v8/blank) + plot_layout(heights=c(0.35,1,0.4,1,0.4,1,0.4,1,0.35))
pdf(file=paste0(output,"Tac1DEAs_col1b.pdf"),height=17, width=8)
(col1b|col2CPM) + plot_layout(widths=c(1,0.2))
dev.off()

#pdf(file=paste0(output,"Tac1DEAs_col1c.pdf"),height=17)
#(v3/blank/v9/blank/v12/blank/v10) + plot_layout(heights=c(1,0.4,1,0.4,1,0.4,1))
#dev.off()
col1c <- (blank/v3/blank/v9/blank/v12/blank/v10/blank) + plot_layout(heights=c(0.35,1,0.4,1,0.4,1,0.4,1,0.35))
pdf(file=paste0(output,"Tac1DEAs_col1c.pdf"),height=17, width=8)
(col1c|col3CPM) + plot_layout(widths=c(1,0.2))
dev.off()



#pdf(file=paste0(output,"Tac1DEAs_col1a.pdf"),height=17)
#(v1/blank/v5/blank/v7/blank/v11) + plot_layout(heights=c(1,0.4,1,0.4,1,0.4,1))
#dev.off()
col1a <- (blank/v1/blank/v5/blank/v7/blank/v11/blank)+ plot_layout(heights=c(0.35,1,0.4,1,0.4,1,0.4,1,0.35))
pdf(file=paste0(output,"Tac1DEAs_col1a.pdf"),height=17, width=8)
(col1a|col1CPM) + plot_layout(widths=c(1,0.2))
dev.off() 

col4a <- (blank/v17/blank/v18/blank/v17/blank/v18/blank) + plot_layout(heights=c(0.35,1,0.4,1,0.4,1,0.4,1,0.35))
#pdf(file=paste0(output,"Tac1DEAs_col4a.pdf"),height=17)
#(v17/blank/v18/blank/blank/blank/blank) + plot_layout(heights=c(0.35,1,0.4,1,0.4,1,0.4,1,0.35))
#dev.off()
pdf(file=paste0(output,"Tac1DEAs_col4a.pdf"),height=17, width=8)
(col4a|col4CPM)+ plot_layout(widths=c(1,0.2))
dev.off()
counts <- mySeurad@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes <- c("Tac1","Adra1a","Greb1","Sox5","Synpr","Pura","Socs2","Mical2","Cck","Phf21b","Fam19a2","Kcnh7","Tiparp","Col25a1","Etv1","Maml2","B3gnt2","Nfib","Moxd1","Rreb1","Syne2")
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
