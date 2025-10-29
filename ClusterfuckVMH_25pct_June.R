#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/VMH_IndependentAnalysis/VMH_Merged_Clusterfuck_JuneFinal_25pct_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")


genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/VMH_MvP_1.5.txt")
eDEG <- read.table("Genelists/VMH_PvU_1.5.txt")
sDEG2 <- read.table("Genelists/VMH_MvU_1.5.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
allSDEGs <- union(Dimorphic, sDEG2)
allSDEGs <- as.matrix(allSDEGs)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)
allSDEGs.filtered <- intersect(allSDEGs, genes.10x)
allSDEGs.filtered <- as.matrix(allSDEGs.filtered)
sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)


AllsDegs <- read.table("topGO/Total_ByRegion/VMH_Genesonly.txt")
AllsDegs <- unlist(AllsDegs)
AllsDegs <- as.matrix(AllsDegs)
AllsDegs.filtered <- intersect(AllsDegs, genes.10x)
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)
totdegs <- nrow(AllsDegs)
detecdegs <- nrow(AllsDegs.filtered)
conc <- cbind(detecdegs, totdegs)
conc.df <- data.frame(conc)
write.csv(conc.df, file=paste0(output,"10xtrapconcordance.csv"))
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)


mySeurat <- NormalizeData(mySeurat, normalization.method="RC", scale.factor=1e6)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


Pct <- DotPlot(mySeurat, features=Dimorphic.filtered)

data <- Pct$data
head(data)
data
data.df <- data.frame(data)
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), mvpprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MvP_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"MvP_clustersperDEG.csv"))


Pct <- DotPlot(mySeurat, features=eDEG.filtered)

data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), pvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"PvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), pvuprop=num/n_clust)
write.csv(dist, file=paste0(output,"PvU_clustersperDEG.csv"))

Pct <- DotPlot(mySeurat, features=sDEG2.filtered)

data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), mvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), mvuprop=num/n_clust)
write.csv(dist, file=paste0(output,"MvU_clustersperDEG.csv"))

Pct <- DotPlot(mySeurat, features=allSDEGs.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), mvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MalevFemale_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), mvuprop=num/n_clust)
write.csv(dist, file=paste0(output,"MalevFemale_clustersperDEG.csv"))

genestotal=data.frame()
for (g in AllsDegs.filtered)
try({

pct <- sum(GetAssayData(mySeurat, slot="data")[g,]>0)/nrow(mySeurat@meta.data)
df <- data.frame(g, pct)
genestotal=rbind(genestotal, df)
})
write.csv(genestotal, file=paste0(output,"MvU_allcellspct.csv"))


Bigplot <- DotPlot(mySeurat, features=AllsDegs.filtered)

data <- Bigplot$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), totalprop=num/n_clust)
head(dist)
allclust <- dist[dist$num == 23,]
allclust
write.csv(allclust, file=paste0(output,"allubiquitousDEGs.csv"))
unique <- dist[dist$num == 1, ]
unique
uniqueclust <- unique$features.plot
uniqueclust.df <- data.frame(uniqueclust)
uniqueclustgenes <- data.df %>% filter(features.plot %in% uniqueclust)
uniqueclustassign <- uniqueclustgenes[uniqueclustgenes$pct.exp > 25,]
uniqueclustassign
write.csv(uniqueclustassign, file=paste0(output,"VMH_clusterspecDEGs_25pct.csv"))
