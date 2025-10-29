library(Seurat)
library(tidyverse)
library(dplyr)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_20pcs_0.1_1.5cutoff"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_20pcs.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
MvPsdegs <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
unlist(MvPsdegs)
MvPsdegs <- as.matrix(MvPsdegs)
MvPsdegs.filtered <- intersect(MvPsdegs, genes.10x)
MvUsdegs <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
unlist(MvUsdegs)
MvUsdegs <- as.matrix(MvUsdegs)
MvUsdegs.filtered <- intersect(MvUsdegs, genes.10x)
PvUsdegs <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
unlist(PvUsdegs)
PvUsdegs <- as.matrix(PvUsdegs)
PvUsdegs.filtered <- intersect(PvUsdegs, genes.10x)
subset.matrixMvP <- mySeurat[MvPsdegs.filtered,]

SFARI <- read.table("Genelists/SFARI.txt")
unlist(SFARI)
SFARI <- as.matrix(SFARI)
SFARI.filtered <- intersect(SFARI, genes.10x)
subset.matrixSFARI <- mySeurat[SFARI.filtered,]

pdf(file=paste0(output,"_SFARIdimorphic.pdf"))
FeaturePlot(subset.matrixSFARI, features=c("nFeature_RNA"), cols=c("light blue", "dark red"), order=TRUE)
dev.off()

pdf(file=paste0(output,"_MvP_binarizedUMAP.pdf"))
FeaturePlot(subset.matrixMvP, features=c("nFeature_RNA"), cols=c("light blue", "dark red"), order=TRUE)
dev.off()
pdf(file=paste0(output,"MvP_binarizedVln.pdf"))
VlnPlot(subset.matrixMvP, features=c("nFeature_RNA"), pt.size=0, ncol=1)
dev.off()
subset.matrixMvU <- mySeurat[MvUsdegs.filtered,]
pdf(file=paste0(output,"_MvP_binarizedUMAP.pdf"))
FeaturePlot(subset.matrixMvP, features=c("nFeature_RNA"), cols=c("light blue", "red"), order=TRUE)
dev.off()
pdf(file=paste0(output,"MvP_binarizedVln.pdf"))
VlnPlot(subset.matrixMvP, features=c("nFeature_RNA"), pt.size=0, ncol=1)
dev.off()
subset.matrixMvU <- mySeurat[MvUsdegs.filtered,]
pdf(file=paste0(output,"MvU_binarizedUMAP.pdf"))
FeaturePlot(subset.matrixMvU, features=c("nFeature_RNA"), cols=c("orange", "green"), order=TRUE)
dev.off()
pdf(file=paste0(output,"MvU_binarizedVln.pdf"))
VlnPlot(subset.matrixMvU, features=c("nFeature_RNA"), pt.size=0, ncol=1)
dev.off()
subset.matrixPvU <- mySeurat[PvUsdegs.filtered,]
pdf(file=paste0(output,"PvU_binarizedUMAP.pdf"))
FeaturePlot(subset.matrixPvU, features=c("nFeature_RNA"), cols=c("yellow", "blue"), order=TRUE)
dev.off()
pdf(file=paste0(output,"PvU_binarizedVln.pdf"))
VlnPlot(subset.matrixPvU, features=c("nFeature_RNA"), pt.size=0, ncol=1)
dev.off()

data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)
avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=MvPsdegs.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in MvPsdegs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "MvP_Clusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "MvP_SDEGS_Per_Cluster.csv"))

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=MvUsdegs.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in MvUsdegs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "MvU_Clusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "MvU_SDEGS_Per_Cluster.csv"))

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=PvUsdegs.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in PvUsdegs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "PvU_Clusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "PvU_SDEGS_Per_Cluster.csv"))

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=SFARI.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in SFARI.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "SFARI_Clusters_Per_Gene.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "SFARI_Genes_Per_Cluster.csv"))