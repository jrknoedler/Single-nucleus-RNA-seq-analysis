#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/PvU_Inhib/PvUinhib_DEGome1.5fold"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
mySeurat <- subset(mySeurat, idents=c("1","3","4","5","6","10","11","12","13","14","17","19","20","21","23","24","25","26","27","29","31","32","33","34","35"))
inhib <- c("1","3","4","5","6","10","11","12","13","14","17","19","20","21","23","24","25","26","27","29","31","32","33","34","35")

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
genelist <- read.table("Genelists/POA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
total_dimorphic = data.frame()
master_list = data.frame()
for (i in inhib){
try({
ident <- i

Male <- subset(x=dsSeurat, subset=Hormone=="Primed")
head(Male[[]])
Male <- subset(x=Male, idents=c(i))
Female <- subset(x=dsSeurat, subset=Hormone=="Unprimed")
Female <- subset(x=Female, idents=c(i))
Maledata <- Male@assays[["RNA"]]@counts
Male.mat <- as.matrix(Maledata)
MaleNum <- ncol(Male.mat)
Femaledata <- Female@assays[["RNA"]]@counts
Female.mat <- as.matrix(Femaledata)
FemaleNum <- ncol(Female.mat)
Combined <- cbind(Male.mat, Female.mat)
ncol(Combined)
cond1 <- "Intact"
cond2 <- "Primed"
sampleTable <- data.frame(condition=factor(c(rep(cond1, MaleNum), rep(cond2, FemaleNum))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
restargeted <- res[filtered.genelist,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(restargeted, file=paste0(output,i,"_sDEGpadjust.csv"))
restargeted <- na.omit(restargeted)
padjcounts <- nrow(restargeted[restargeted$padj < 0.05,])
df <- data.frame(i, padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- restargeted[restargeted$padj < 0.05,]
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(output,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(output,"totalsigDEGs.csv"))
genes <- rownames(master_list)
genes <- unlist(genes)
genes <- unique(genes)
write.csv(genes, file=paste0(output, "uniquesigDEGs.csv"))
