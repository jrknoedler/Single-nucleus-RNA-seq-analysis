#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Marchredo_1cpm_"

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
mySeurat <- NormalizeData(mySeurat, normalization.method="RC", scale.factor=1e6)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=Dimorphic.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in Dimorphic.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "MvPClusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 1)
write.csv(counts, file=paste0(output, "MvPSDEGS_Per_Cluster.csv"))

avg2 <- merged %>% group_by(clusters) %>% summarize_at(.vars=eDEG.filtered, .funs=c("mean"))
avg2
avg2 %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total2= data.frame()
for (g in eDEG.filtered)
try({
results2 <- paste0(g,"results")
results2 = NULL
counts2 <- nrow(avg2[avg2[[g]] > 1,])
df2 <- data.frame(g, counts2)
df_total2 = rbind(df_total2, df2)
})

write.csv(df_total2, file=paste0(output, "PvUClusters_Per_SDEG.csv"))
avg2 <- data.matrix(avg2)
counts2 <- rowSums(avg2 > 1)
write.csv(counts2, file=paste0(output, "PvUSDEGS_Per_Cluster.csv"))

avg3 <- merged %>% group_by(clusters) %>% summarize_at(.vars=sDEG2.filtered, .funs=c("mean"))
avg3
avg3 %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total3= data.frame()
for (g in sDEG2.filtered)
try({
results3 <- paste0(g,"results")
results3 = NULL
counts3 <- nrow(avg3[avg3[[g]] > 1 ,])
df3 <- data.frame(g, counts3)
df_total3 = rbind(df_total3, df3)
})

write.csv(df_total3, file=paste0(output, "MvUClusters_Per_SDEG.csv"))
avg3 <- data.matrix(avg3)
counts3 <- rowSums(avg3 > 1)
write.csv(counts3, file=paste0(output, "MvUSDEGS_Per_Cluster.csv"))