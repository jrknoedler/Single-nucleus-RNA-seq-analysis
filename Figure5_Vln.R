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

##Load VMH and compile violin data
output <- "Seurat/VMH_IndependentAnalysis/VMH_Figs5heatmapsfinal_finalrecluster"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
mySeurat <-NormalizeData(mySeurat)

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)
mylevels <- rev(mylevels)

dat <- VlnPlot(mySeurat, features=c("Cckar"))
data <- dat$data
data$val <- "0"



p1 <- ggplot(data, aes(x=Cckar, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

pdf(file=paste0(output,"Cckarhorizonatltest.pdf"), height=30, width=5)
p1
dev.off()

dat <- VlnPlot(mySeurat, features=c("Tnfaip8l3"))
data <- dat$data
data$val <- "0"



p2 <- ggplot(data, aes(x=Tnfaip8l3, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(mySeurat, features=c("Mme"))
data <- dat$data
data$val <- "0"



p3 <- ggplot(data, aes(x=Mme, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(mySeurat, features=c("Cgnl1"))
data <- dat$data
data$val <- "0"



p4 <- ggplot(data, aes(x=Cgnl1, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(mySeurat, features=c("Popdc3"))
data <- dat$data
data$val <- "0"



p5 <- ggplot(data, aes(x=Popdc3, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(mySeurat, features=c("Trim36"))
data <- dat$data
data$val <- "0"



p6 <- ggplot(data, aes(x=Trim36, y=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

pdf(file=paste0(output,"VMHFig5Vln.pdf"), width=40, height=70)
(p1+p2+p3+p4+p5+p6) + plot_layout(nrow=1) & theme(axis.text=element_blank())
dev.off()




##Load BNST and compile violin data
BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
BNST <- NormalizeData(BNST)

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
BNSTlevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)
BNSTlevels <- rev(BNSTlevels)
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)

dat <- VlnPlot(BNST, features=c("Rbm20"))
data <- dat$data
data$val <- "0"



bp1 <- ggplot(data, aes(x=Rbm20, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Plekhg1"))
data <- dat$data
data$val <- "0"


bp2 <- ggplot(data, aes(x=Plekhg1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Prok2"))
data <- dat$data
data$val <- "0"


bp3 <- ggplot(data, aes(x=Prok2, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Cck"))
data <- dat$data
data$val <- "0"


bp4 <- ggplot(data, aes(x=Cck, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Phf21b"))
data <- dat$data
data$val <- "0"


bp5 <- ggplot(data, aes(x=Phf21b, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Slco2a1"))
data <- dat$data
data$val <- "0"


bp6 <- ggplot(data, aes(x=Slco2a1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Tac1"))
data <- dat$data
data$val <- "0"


bp7 <- ggplot(data, aes(x=Tac1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Asic4"))
data <- dat$data
data$val <- "0"


bp8 <- ggplot(data, aes(x=Asic4, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Kank1"))
data <- dat$data
data$val <- "0"


bp9 <- ggplot(data, aes(x=Kank1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Sst"))
data <- dat$data
data$val <- "0"


bp10 <- ggplot(data, aes(x=Sst, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Abca1"))
data <- dat$data
data$val <- "0"


bp11 <- ggplot(data, aes(x=Abca1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Amot"))
data <- dat$data
data$val <- "0"


bp12 <- ggplot(data, aes(x=Amot, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Fam107b"))
data <- dat$data
data$val <- "0"


bp13 <- ggplot(data, aes(x=Fam107b, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Prkcd"))
data <- dat$data
data$val <- "0"


bp14 <- ggplot(data, aes(x=Prkcd, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Sertm1"))
data <- dat$data
data$val <- "0"


bp15 <- ggplot(data, aes(x=Sertm1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Tns1"))
data <- dat$data
data$val <- "0"


bp16 <- ggplot(data, aes(x=Tns1, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Sytl4"))
data <- dat$data
data$val <- "0"


bp17 <- ggplot(data, aes(x=Sytl4, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))


dat <- VlnPlot(BNST, features=c("Npr3"))
data <- dat$data
data$val <- "0"


bp18 <- ggplot(data, aes(x=Npr3, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(BNST, features=c("Plscr4"))
data <- dat$data
data$val <- "0"


bp19 <- ggplot(data, aes(x=Plscr4, y=factor(ident, levels=BNSTlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
blank <- ggplot() +theme_void()

##test an intermediate plot
BNST <- (bp1|bp2|bp3|bp4|bp5|bp6|bp7|bp8|bp9|bp10|bp11|bp12|bp13|bp14|bp15|bp16|bp17|bp18|bp19) + plot_layout(nrow=1) + theme(axis.text=element_blank())
pdf(file=paste0(output,"BNST_test.pdf"), width=50, height=30)
BNST & theme(axis.text=element_blank())
dev.off()
VMH <- (p1|p2|p3|p4|p5|p6) 
pdf(file=paste0(output,"BNST_VMHvln.pdf"), width=100, height=40)
BNST+VMH
dev.off()

MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

MeA <- NormalizeData(MeA)

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
MeAlevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34)
MeAlevels <- rev(MeAlevels)
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)

dat <- VlnPlot(MeA, features=c("Sst"))
data <- dat$data
data$val <- "0"



mp1 <- ggplot(data, aes(x=Sst, y=factor(ident, levels=MeAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(MeA, features=c("Map3k15"))
data <- dat$data
data$val <- "0"



mp2 <- ggplot(data, aes(x=Map3k15, y=factor(ident, levels=MeAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(MeA, features=c("Tspan18"))
data <- dat$data
data$val <- "0"



mp3 <- ggplot(data, aes(x=Tspan18, y=factor(ident, levels=MeAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(MeA, features=c("Vgll3"))
data <- dat$data
data$val <- "0"



mp4 <- ggplot(data, aes(x=Vgll3, y=factor(ident, levels=MeAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(MeA, features=c("Efemp1"))
data <- dat$data
data$val <- "0"



mp5 <- ggplot(data, aes(x=Efemp1, y=factor(ident, levels=MeAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

pdf(file=paste0(output,"MeAFig5Vln.pdf"), width=30, height=50)
(mp1+mp2+mp3+mp4+mp5) + plot_layout(nrow=1) & theme(axis.text=element_blank())
dev.off()


POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")

POA <- NormalizeData(POA)

POAlevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)
POAlevels <- rev(POAlevels)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)

dat <- VlnPlot(POA, features=c("Efemp1"))
data <- dat$data
data$val <- "0"



pp1 <- ggplot(data, aes(x=Efemp1, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))


dat <- VlnPlot(POA, features=c("Mob3b"))
data <- dat$data
data$val <- "0"



pp2 <- ggplot(data, aes(x=Mob3b, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Moxd1"))
data <- dat$data
data$val <- "0"



pp3 <- ggplot(data, aes(x=Moxd1, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("St3gal1"))
data <- dat$data
data$val <- "0"



pp4 <- ggplot(data, aes(x=St3gal1, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Thbs4"))
data <- dat$data
data$val <- "0"



pp5 <- ggplot(data, aes(x=Thbs4, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Mrc2"))
data <- dat$data
data$val <- "0"



pp6 <- ggplot(data, aes(x=Mrc2, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Ttn"))
data <- dat$data
data$val <- "0"



pp7 <- ggplot(data, aes(x=Ttn, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Neb"))
data <- dat$data
data$val <- "0"



pp8 <- ggplot(data, aes(x=Neb, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Kiss1"))
data <- dat$data
data$val <- "0"



pp9 <- ggplot(data, aes(x=Kiss1, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Kl"))
data <- dat$data
data$val <- "0"



pp10 <- ggplot(data, aes(x=Kl, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Slc17a8"))
data <- dat$data
data$val <- "0"



pp11 <- ggplot(data, aes(x=Slc17a8, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Abtb2"))
data <- dat$data
data$val <- "0"



pp12 <- ggplot(data, aes(x=Abtb2, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

dat <- VlnPlot(POA, features=c("Penk"))
data <- dat$data
data$val <- "0"



pp13 <- ggplot(data, aes(x=Penk, y=factor(ident, levels=POAlevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))

pdf(file=paste0(output,"POAFig5Vln.pdf"), width=30, height=50)
(pp1+pp2+pp3+pp4+pp5+pp6+pp7+pp8+pp9+pp10+pp11+pp12+pp13) + plot_layout(nrow=1) & theme(axis.text=element_blank())
dev.off()
