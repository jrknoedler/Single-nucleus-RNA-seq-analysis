#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/Reclustered_DEGome/Default/MvP/MvP_DEGome1.5fold_reclusteredfinal"
mvuoutput <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/Reclustered_DEGome/Default/MvU/MvU_DEGome1.5fold_reclusteredfinal"
pvuoutput <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/Reclustered_DEGome/Default/PvU/PvU_DEGome1.5fold_reclusteredfinal"
baseoutput <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/Reclustered_DEGome/Default/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)



datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone


dsSeurat <- AddMetaData(dsSeurat, metadata=idents, col.name="Ident")
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
Idents(dsSeurat) <- "Ident"
head(dsSeurat[[]])
dsSeurat <- subset(dsSeurat, nCount_RNA > 7500)
genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvP <- intersect(genelist, genes.10x)
filtered.MvP
dsSeurat$celltype.Hormone <- paste(Idents(dsSeurat), dsSeurat$Hormone, sep="_")
dsSeurat$celltype <- Idents(dsSeurat)
Idents(dsSeurat) <- "celltype.Hormone"
total_dimorphic = data.frame()
master_list = data.frame()
DefaultAssay(mySeurat) <- "RNA"
dsSeurat <- NormalizeData(dsSeurat)
for (i in 0:35){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(dsSeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.MvP, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(output,i,"_MASTDEGOME.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
padjcounts <- nrow(sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,])
df <- data.frame(i,padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(output,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(output,"totalsigDEGs.csv"))
genes <- rownames(master_list)
genes <-  master_list$gene
genes <- unique(genes)
write.csv(genes, file=paste0(output, "uniquesigDEGs.csv"))


genelist <- read.table("Genelists/BNST_MvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvU <- intersect(genelist, genes.10x)
filtered.MvU
total_dimorphic = data.frame()
master_list = data.frame()
for (i in 0:35){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Unprimed")
sex.dimorphism <- FindMarkers(dsSeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE,features=filtered.MvU, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(mvuoutput,i,"_MASTDEGOME.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
padjcounts <- nrow(sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,])
df <- data.frame(i,padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(mvuoutput,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(mvuoutput,"totalsigDEGs.csv"))
genes <- master_list$gene
genes <- unlist(genes)
genes <- unique(genes)
write.csv(genes, file=paste0(mvuoutput, "uniquesigDEGs.csv"))


genelist <- read.table("Genelists/BNST_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.PvU <- intersect(genelist, genes.10x)
filtered.PvU
total_dimorphic = data.frame()
master_list = data.frame()
for (i in 0:35){
try({
ident1 <- paste0(i,"_Primed")
ident2 <- paste0(i,"_Unprimed")
sex.dimorphism <- FindMarkers(dsSeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.PvU, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(pvuoutput,i,"_MASTDEGOME.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
padjcounts <- nrow(sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,])
df <- data.frame(i,padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(pvuoutput,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(pvuoutput,"totalsigDEGs.csv"))
genes <- master_list$gene
genes <- unlist(genes)
genes <- unique(genes)
write.csv(genes, file=paste0(pvuoutput, "uniquesigDEGs.csv"))


Idents(dsSeurat) <- "Hormone"
MvP <- FindMarkers(dsSeurat, assay="RNA", ident.1 = "Intact", ident.2="Primed", min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.MvP, verbose=TRUE)
write.csv(MvP, file=paste0(baseoutput,"MvP_Metadegs.csv"))
MvU <- FindMarkers(dsSeurat, assay="RNA", ident.1 = "Intact", ident.2="Unprimed", min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.MvU, verbose=TRUE)
write.csv(MvU, file=paste0(baseoutput,"MvU_Metadegs.csv"))
PvU <- FindMarkers(dsSeurat, assay="RNA", ident.1 = "Primed", ident.2="Unprimed", min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.PvU, verbose=TRUE)
write.csv(PvU, file=paste0(baseoutput,"PvU_Metadegs.csv"))