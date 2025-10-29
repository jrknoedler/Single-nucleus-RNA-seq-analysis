#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Marchredo_25pct_"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
eDEG <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
sDEG2 <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


Pct <- DotPlot(mySeurat, features=Dimorphic.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 5), mvpprop=num/n_cells)
pct

dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp > 5), mvpprop=num/n_clust)



avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=Dimorphic.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")


logtest <- avg[,"Pappa"]
head(logtest)
datatest <- data[data$features.plot=="Pappa",]
head(datatest)
combined <- cbind(logtest,datatest)
combined
combinedsub <- combined[combined$pct.exp > 5 & combined$Pappa > 0.1,]
combinedsub

df_total= data.frame()
for (g in Dimorphic.filtered)

{
logcounts <- avg[,g]
colnames(logcounts) <- c("gene")
props  <- data[data$features.plot==g,]
combined <- cbind(logcounts,props)
saveRDS(combined, file=paste0(output,"gene.rds"))
count <- nrow(combined[combined$pct.exp >5 & combined$gene > 0.1,])
mCTcount <- data.frame(g, count)
df_total=rbind(df_total, mCTcount)
}
write.csv(df_total, file=paste0(output,"IntersectCountMvP.csv"))

avg2 <- merged %>% group_by(clusters) %>% summarize_at(.vars=eDEG.filtered, .funs=c("mean"))
avg2
avg2 %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total2= data.frame()
for (g in eDEG.filtered)
try({
results2 <- paste0(g,"results")
results2 = NULL
counts2 <- nrow(avg2[avg2[[g]] > 0.1,])
df2 <- data.frame(g, counts2)
df_total2 = rbind(df_total2, df2)
})



Pct <- DotPlot(mySeurat, features=eDEG.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 5), pvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"PvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >5), pvuprop=num/n_clust)
write.csv(dist, file=paste0(output,"PvU_clustersperDEG.csv"))

df_total= data.frame()
for (g in eDEG.filtered)

{
logcounts <- avg2[,g]
colnames(logcounts) <- c("gene")
props  <- data[data$features.plot==g,]
combined <- cbind(logcounts,props)
saveRDS(combined, file=paste0(output,"gene.rds"))
count <- nrow(combined[combined$pct.exp >5 & combined$gene > 0.1,])
mCTcount <- data.frame(g, count)
df_total=rbind(df_total, mCTcount)
}
write.csv(df_total, file=paste0(output,"IntersectCountPvU.csv"))



avg3 <- merged %>% group_by(clusters) %>% summarize_at(.vars=sDEG2.filtered, .funs=c("mean"))
avg3
avg3 %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total3= data.frame()
for (g in sDEG2.filtered)
try({
results3 <- paste0(g,"results")
results3 = NULL
counts3 <- nrow(avg3[avg3[[g]] > 0.1 ,])
df3 <- data.frame(g, counts3)
df_total3 = rbind(df_total3, df3)
})


Pct <- DotPlot(mySeurat, features=sDEG2.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 5), mvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >5), mvuprop=num/n_clust)
write.csv(dist, file=paste0(output,"MvU_clustersperDEG.csv"))

df_total= data.frame()
for (g in sDEG2.filtered)

{
logcounts <- avg3[,g]
colnames(logcounts) <- c("gene")
props  <- data[data$features.plot==g,]
combined <- cbind(logcounts,props)
saveRDS(combined, file=paste0(output,"gene.rds"))
count <- nrow(combined[combined$pct.exp >5 & combined$gene > 0.1,])
mCTcount <- data.frame(g, count)
df_total=rbind(df_total, mCTcount)
}
write.csv(df_total, file=paste0(output,"IntersectCountMvU.csv"))