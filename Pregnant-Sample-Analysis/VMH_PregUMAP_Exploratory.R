#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(scales)
output <- "X:/Seurat/VMH_IndependentAnalysis/VMH_Merged_PregDEGProject"
memory.limit(64000)
mySeurat <- readRDS("X:/Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
PregvMale <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvMale_all1.5.txt")
PregvEst <- read.table("X:/Genelists/PregnancyDEGs/VMH_PvFr_All1.5.txt")
PregvDiest <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvFnr_All1.5.txt")
unlist(PregvEst)
unlist(PregvMale)
unlist(PregvDiest)
PregvEst.m <- as.matrix(PregvEst)
PregvEst.filtered <- intersect(PregvEst.m, genes.10x)
PregvEst.absent <- setdiff(PregvEst.m,genes.10x)
PregvEst.filtered <- as.matrix(PregvEst.filtered)
PregvMale.m <- as.matrix(PregvMale)
PregvMale.filtered <- intersect(PregvMale.m, genes.10x)
PregvMale.absent <- setdiff(PregvMale.m, genes.10x)
PregvMale.filtered <- as.matrix(PregvMale.filtered)
PregvDiest.m <- as.matrix(PregvDiest)
PregvDiest.filtered <- intersect(PregvDiest.m, genes.10x)
PregvDiest.absent <- setdiff(PregvDiest.m, genes.10x)
PregvDiest.filtered <- as.matrix(PregvDiest.filtered)
allDEGs1 <- union(PregvMale, PregvEst)
allDEGs <- union(allDEGs1, PregvDiest)
allDEGs <- as.matrix(allDEGs)
allDEGs.filtered <- intersect(allDEGs, genes.10x)

PregAbsent1 <- union(PregvDiest.absent,PregvEst.absent)
PregAbsent <- union(PregAbsent1, PregvMale.absent)
write.csv(PregAbsent,file=paste0(output,"_PregDEGsnotfound.csv"))
Pct <- DotPlot(mySeurat, features=PregvEst.filtered)

data <- Pct$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(MvPreg=sum(pct.exp > 10), mvpregprop=MvPreg/n())
head(pct)
pct <- column_to_rownames(pct, 'id')
head(pct)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"_PregvEst_clustersperDEG.csv"))


Pct <- DotPlot(mySeurat, features=PregvMale.filtered)

data <- Pct$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(MvPreg=sum(pct.exp > 10), mvpregprop=MvPreg/n())
head(pct)
pct <- column_to_rownames(pct, 'id')
head(pct)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"_PregvMale_clustersperDEG.csv"))


Pct <- DotPlot(mySeurat, features=PregvDiest.filtered)

data <- Pct$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(MvPreg=sum(pct.exp > 10), mvpregprop=MvPreg/n())
head(pct)
pct <- column_to_rownames(pct, 'id')
head(pct)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"_PregvDiest_clustersperDEG.csv"))




allDEGs1 <- union(PregvMale, PregvDiest)
allDEGs <- union(allDEGs1,PregvEst)
allDEGs <- as.matrix(allDEGs)
allDEGs.filtered <- intersect(allDEGs, genes.10x)
allDEGs.filtered <- as.matrix(allDEGs.filtered)

Up_MvPreg <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvMale_DownPreg.txt", header=FALSE)
Up_MvPreg <- unlist(Up_MvPreg)
Up_MvPreg.m <- as.matrix(Up_MvPreg)
Up_MvPreg.filtered <- intersect(Up_MvPreg.m, genes.10x)

Up_PregvM <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvMale_UpPreg.txt", header=FALSE)
Up_PregvM <- unlist(Up_PregvM)
Up_PregvM.m <- as.matrix(Up_PregvM)
Up_PregvM.filtered <- intersect(Up_PregvM.m, genes.10x)

Up_PregvEst <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvFr_UpPreg.txt", header=FALSE)
Up_PregvEst <- unlist(Up_PregvEst)
Up_PregvEst.m <- as.matrix(Up_PregvEst)
Up_PregvEst.filtered <- intersect(Up_PregvEst.m, genes.10x)

Up_EstvPreg <- read.table("X:/Genelists/PregnancyDEGs/VMH_PregvFr_DownPreg.txt", header=FALSE)
Up_EstvPreg <- unlist(Up_EstvPreg)
Up_EstvPreg.m <- as.matrix(Up_EstvPreg)
Up_EstvPreg.filtered <- intersect(Up_EstvPreg.m, genes.10x)

Up_PregvDiest <- read.table("X:Genelists/PregnancyDEGs/VMH_PregvFnr_UpPreg.txt", header=FALSE)
Up_PregvDiest <- unlist(Up_PregvDiest)
Up_PregvDiest.m <- as.matrix(Up_PregvDiest)
Up_PregvDiest.filtered <- intersect(Up_PregvDiest.m, genes.10x)

Up_DiestvPreg <- read.table("X:Genelists/PregnancyDEGs/VMH_PregvFnr_DownPreg.txt", header=FALSE)
Up_DiestvPreg <- unlist(Up_DiestvPreg)
Up_DiestvPreg.m <- as.matrix(Up_DiestvPreg)
Up_DiestvPreg.filtered <- intersect(Up_DiestvPreg.m, genes.10x)

PregUpAll <- union(Up_PregvDiest,Up_PregvEst)
PregUpUni <- intersect(Up_PregvDiest,Up_PregvEst)
PregUpUni <- unlist(PregUpUni)
write.csv(PregUpUni, file=paste0(output,"AllCommonUpPreg.csv"))
PregUpUni.m <- as.matrix(PregUpUni)
PregUpUni.filtered <- intersect(genes.10x, PregUpUni.m)
write.csv(PregUpUni.filtered, file=paste0(output,"_PregUpEst_AND_Diest.csv"))

PregDownAll <- union(Up_DiestvPreg,Up_EstvPreg)
PregDownUni <- intersect(Up_DiestvPreg,Up_EstvPreg)
write.csv(PregDownUni, file=paste0(output,"AllCommonDownPreg.csv"))
PregDownUni <- unlist(PregDownUni)
PregDownUni.m <- as.matrix(PregDownUni)
PregDownUni.filtered <- intersect(genes.10x, PregDownUni.m)
write.csv(PregDownUni.filtered, file=paste0(output,"_PregDownEst_AND_Diest.csv"))

Femalecompup <- unlist(PregUpAll)
Femalecompup.m <- as.matrix(Femalecompup)
Femalecompup.filtered <- intersect(Femalecompup.m, genes.10x)
PregupFin <- union(PregUpAll,Up_PregvM)
PregupFin <- unlist(PregupFin)
PregupFin.m <- as.matrix(PregupFin)
PregupFin.filtered <- intersect(PregupFin.m, genes.10x)

PregDownAll <- union(Up_DiestvPreg,Up_EstvPreg)
Femalecompdown <- unlist(PregDownAll)
Femalecompdown.m <- as.matrix(Femalecompdown)
Femalecompdown.filtered <- intersect(Femalecompdown.m, genes.10x)
PregDownFin <- union(PregDownAll,Up_MvPreg)
PregDownFin <- unlist(PregDownFin)
PregDownFin.m <- as.matrix(PregDownFin)
PregDownFin.filtered <- intersect(PregDownFin.m, genes.10x)


data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)

PregupNamed=PregupFin.filtered
PregdownNamed=PregDownFin.filtered
FemaleCompUpnamed=Femalecompup.filtered
FemaleCompDownnamed=Femalecompdown.filtered
PregUpUniNamed=PregUpUni.filtered
PregDownUniNamed=PregDownUni.filtered
PregupNamed=list("Pregup"=as.character(PregupNamed))
PregdownNamed=list("Pregdown"=as.character(PregdownNamed))
FemalePregupnamed=list("F_PregUp"=as.character(FemaleCompUpnamed))
FemalePregdownnamed=list("F_PregDown"=as.character(FemaleCompDownnamed))
FemalePregUpUni=list("F_PregSignature"=as.character(PregUpUniNamed))
FemalePregDownUni=list("F_PregdownSignature"=as.character(PregDownUniNamed))

Pctuni <- DotPlot(mySeurat, features=PregUpUni.filtered)
data <- Pctuni$data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
Pctuni <- data.df %>% group_by(id) %>% summarize(PregDown=sum(pct.exp>10), PregDownProp=PregDown/n())

clusters$PregDownProp <- Pctpd$PregDownProp[match(clusters$clusters,Pctpd$id)]
clustersPregDown <- clusters["PregDownProp"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregDown, col.name="IntersectPregUpDEGs")

pdf(file=paste0(output,"_Intersect_PregUpvF.pdf"))
FeaturePlot(mySeurat, features=c("IntersectPregUpDEGs"))
dev.off()

Pctuni <- DotPlot(mySeurat, features=PregDownUni.filtered)
data <- Pctuni$data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
Pctuni <- data.df %>% group_by(id) %>% summarize(PregDown=sum(pct.exp>10), PregDownProp=PregDown/n())

clusters$PregDownProp <- Pctpd$PregDownProp[match(clusters$clusters,Pctpd$id)]
clustersPregDown <- clusters["PregDownProp"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregDown, col.name="IntersectPregUpDEGs")

pdf(file=paste0(output,"_Intersect_PregDownvF.pdf"))
FeaturePlot(mySeurat, features=c("IntersectPregUpDEGs"))
dev.off()


pdf(file=paste0(output,"_PregUpUnionDotplot.pdf"), height=14, width=70)
DotPlot(mySeurat, features=PregUpUni.filtered)
dev.off()

pdf(file=paste0(output,"_PregDownUnionDotplot.pdf"), height=14, width=70)
DotPlot(mySeurat, features=PregDownUni.filtered)
dev.off()


Pctuni <- column_to_rownames(Pctuni, 'id')
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"_IntersectPregUpvFemale_clustersperDEG.csv"))




mySeurat <- AddModuleScore(mySeurat, features=c(PregupNamed,PregdownNamed,FemalePregupnamed,FemalePregdownnamed,FemalePregUpUni,FemalePregDownUni),pool=NULL,seed=1, name=c("PregUpMod","PregDownMod","F_PregupMod","F_PregDownMod","F_UnionPregUp","F_UnionPregDown"))
pdf(file=paste0(output,"_PregMods.pdf"), width=30)
FeaturePlot(mySeurat, features=c("PregUpMod1","PregDownMod2","F_PregupMod3","F_PregDownMod4","F_UnionPregUp5","F_UnionPregDown6"))
dev.off()
pdf(file=paste0(output,"_PregMod_Vln.pdf"), width=28, height=20)
VlnPlot(mySeurat, features=c("PregUpMod1","PregDownMod2","F_PregupMod3","F_PregDownMod4","F_UnionPregUp5","F_UnionPregDown6"), pt.size=0, ncol=1)
dev.off()

Pctpu <- DotPlot(mySeurat, features=PregupFin.filtered)
data <- Pctpu$data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
Pctpu <- data.df %>% group_by(id) %>% summarize(PregUp=sum(pct.exp>10), PregUpProp=PregUp/n())

clusters$PregUpProp <- Pctpu$PregUpProp[match(clusters$clusters,Pctpu$id)]
clustersPregUp <- clusters["PregUpProp"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregUp, col.name="PregUpDEGs")

pdf(file=paste0(output,"_PregUpUmap.pdf"))
FeaturePlot(mySeurat, features="PregUpDEGs")
dev.off()

Pctpd <- DotPlot(mySeurat, features=PregDownFin.filtered)
data <- Pctpd$data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
Pctpd <- data.df %>% group_by(id) %>% summarize(PregDown=sum(pct.exp>10), PregDownProp=PregDown/n())

clusters$PregDownProp <- Pctpd$PregDownProp[match(clusters$clusters,Pctpd$id)]
clustersPregDown <- clusters["PregDownProp"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregDown, col.name="PregDownDEGs")

pdf(file=paste0(output,"PregDownUMAP.pdf"))
FeaturePlot(mySeurat, features="PregDownDEGs")
dev.off()


Pcta <- DotPlot(mySeurat, features=allDEGs.filtered)
data <- Pcta$data

data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pcta <- data.df %>% group_by(id) %>% summarize(AllDEGs=sum(pct.exp > 10), DEGprop=AllDEGs/n())



clusters$DEGprop <- pcta$DEGprop[match(clusters$clusters,pcta$id)]
clustersall <- clusters["DEGprop"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersall, col.name="AllDEGs")

pdf(file=paste0(output,"_AllDEGs.pdf"))
FeaturePlot(mySeurat, features="AllDEGs")
dev.off()

Pct1 <- DotPlot(mySeurat, features=Up_MvPreg.filtered)

data <- Pct1$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct1 <- data.df %>% group_by(id) %>% summarize(MvPreg=sum(pct.exp > 10), mvpregprop=MvPreg/n())
head(pct1)
pct1 <- column_to_rownames(pct1, 'id')
head(pct1)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_MvPregclustersperDEG.csv"))



Pct2 <- DotPlot(mySeurat, features=Up_PregvM.filtered)

data <- Pct2$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct2 <- data.df %>% group_by(id) %>% summarize(PregvM=sum(pct.exp > 10), Pregvmprop=PregvM/n())
head(pct2)
pct2 <- column_to_rownames(pct2, 'id')
head(pct2)
write.csv(pct2, file=paste0(output,"FroverM_sDEGspercluster.csv"))

dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_MvPregclustersperDEG.csv"))


comb <- cbind(pct1, pct2)
head(comb)
PregvMrat <- comb %>% summarize(MPregRatio=Pregvmprop/mvpregprop)


PregvMrat <- scale(PregvMrat)
PregvMrat <- data.frame(PregvMrat)

PregvMrat <- rownames_to_column(PregvMrat, var="id")
head(PregvMrat)


clusters$PregMRatio <- PregvMrat$MPregRatio[match(clusters$clusters,PregvMrat$id)]
clustersMPreg <- clusters["PregMRatio"]
head(clustersMPreg)
mySeurat <- AddMetaData(mySeurat, metadata=clustersMPreg, col.name="MvPregRatio")

#pdf(file=paste0(output,"VMH_MvPregProjecttest.pdf"))
#FeaturePlot(mySeurat, features="MvPregRatio")
#dev.off()

p <- FeaturePlot(mySeurat, features="MvPregRatio")
pdata3 <- p$data
pdf(file=paste0(output,"_VMH_MvPregProjecttest.pdf"))
ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=MvPregRatio)) + cowplot::theme_cowplot() + geom_point(size=0.3)+scale_colour_gradient2(low="blue", high="darkorange1", mid="light gray", guide="colorbar")
dev.off()




Pct3 <- DotPlot(mySeurat, features=Up_PregvEst.filtered)

data <- Pct3$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct3 <- data.df %>% group_by(id) %>% summarize(Pregvest=sum(pct.exp > 10), Pregvestprop=Pregvest/n())
head(pct3)
pct3 <- column_to_rownames(pct3, 'id')
head(pct3)

dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), Pregvestprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_PregvEstclustersperDEG.csv"))


Pct4 <- DotPlot(mySeurat, features=Up_EstvPreg.filtered)

data <- Pct4$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct4 <- data.df %>% group_by(id) %>% summarize(EstvPreg=sum(pct.exp > 10), EstvPregprop=EstvPreg/n())
head(pct4)
pct4 <- column_to_rownames(pct4, 'id')
head(pct4)

dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), Estvpregprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_EstvPregclustersperDEG.csv"))


comb2 <- cbind(pct3, pct4)
comb2
Pestrat <- comb2 %>% summarize(PregEstRatio=Pregvestprop/EstvPregprop)
Pestrat

Pestrat <- scale(Pestrat)
Pestrat <- data.frame(Pestrat)


Pestrat <- rownames_to_column(Pestrat, var="id")
head(Pestrat)


clusters$PregEstRatio <- Pestrat$PregEstRatio[match(clusters$clusters,Pestrat$id)]
clustersPregEst <- clusters["PregEstRatio"]
head(clustersPregEst)
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregEst, col.name="PregEstRatio")

#(file=paste0(output,"VMH_PregvEstProjecttest.pdf"))
#FeaturePlot(mySeurat, features="PregEstRatio")
#dev.off()

p <- FeaturePlot(mySeurat, features="PregEstRatio")
pdata3 <- p$data
pdf(file=paste0(output,"_VMH_PregvEstProjecttest.pdf"))
ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=PregEstRatio)) + cowplot::theme_cowplot() + geom_point(size=0.3)+scale_colour_gradient2(low="deeppink1", high="darkorange1", mid="light gray", guide="colorbar")
dev.off()




Pct5 <- DotPlot(mySeurat, features=Up_PregvDiest.filtered)

data <- Pct5$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct5 <- data.df %>% group_by(id) %>% summarize(PregvDiest=sum(pct.exp > 10), PregvDiestprop=PregvDiest/n())
head(pct5)
pct5 <- column_to_rownames(pct5, 'id')
head(pct5)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_PregvDiestclustersperDEG.csv"))



Pct6 <- DotPlot(mySeurat, features=Up_DiestvPreg.filtered)

data <- Pct6$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct6 <- data.df %>% group_by(id) %>% summarize(DiestvPreg=sum(pct.exp > 10), DiestvPregprop=DiestvPreg/n())
head(pct6)
pct6 <- column_to_rownames(pct6, 'id')
head(pct6)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_DiestvPregclustersperDEG.csv"))



comb <- cbind(pct5, pct6)
head(comb)
PregDierat <- comb %>% summarize(PregdiestRatio=PregvDiestprop/DiestvPregprop)

PregDierat <- scale(PregDierat)
PregDierat <- data.frame(PregDierat)



PregDierat <- rownames_to_column(PregDierat, var="id")
head(PregDierat)


clusters$PregdiestRatio <- PregDierat$PregdiestRatio[match(clusters$clusters,PregDierat$id)]
clustersPregDiest <- clusters["PregdiestRatio"]
head(clustersPregDiest)
mySeurat <- AddMetaData(mySeurat, metadata=clustersPregDiest , col.name="PregDiestratio")

#pdf(file=paste0(output, "_VMH_PregvDiestRatio.pdf"))
#FeaturePlot(mySeurat, features="PregDiestratio", cols=c("chartreuse3","orange"), blend.threshold = 1)
#ev.off()

p <- FeaturePlot(mySeurat, features=c("PregDiestratio"))
pdata3 <- p$data
pdf(file=paste0(output,"_VMH_PregvDiestRatio.pdf"))
ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=PregDiestratio)) + cowplot::theme_cowplot() + geom_point(size=0.3)+scale_colour_gradient2(low="chartreuse3", high="darkorange1", mid="light gray", guide="colorbar")
dev.off()