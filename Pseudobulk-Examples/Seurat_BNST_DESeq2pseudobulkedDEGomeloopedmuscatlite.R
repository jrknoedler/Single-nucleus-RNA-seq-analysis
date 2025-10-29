#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/Unbiased_DEGome/MvP/MvP_UnbiasedDEGome1.5fold_reclustered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(DESeq2)
library(Matrix)
library(matrixStats)



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
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
dsSeurat$celltype.Hormone <- paste(Idents(dsSeurat), dsSeurat$Hormone, sep="_")
dsSeurat$celltype <- Idents(dsSeurat)
Idents(dsSeurat) <- "celltype"
total_dimorphic = data.frame()
master_list = data.frame()
DefaultAssay(mySeurat) <- "RNA"
dsSeurat <- NormalizeData(dsSeurat)

total_dimorphic = data.frame()
master_list = data.frame()

for (i in 0:35){
try({
ident <- i

Male <- subset(x=dsSeurat, subset=Hormone=="Intact")
Male <- subset(x=Male, idents=c(i))
#Male <- subset(Male, subset=nCount_RNA > 7500 & nCount_RNA < 40000)
Female <- subset(x=dsSeurat, subset=Hormone=="Primed")
Female <- subset(x=Female, idents=c(i))
#Female <- subset(Female, subset=nCount_RNA > 7500 & nCount_RNA < 40000)
Maledata <- Male@assays[["RNA"]]@counts
dim(Maledata)
Maledata.df <- as.data.frame(t(Maledata))
dim(Maledata.df)
subgroup <- rep(1:3, ceiling(nrow(Maledata.df)/3))
subgroup
dim(subgroup)
Maledata.df$pseudorep <- 0 
size <- nrow(Maledata.df)
size
Maledata.df[,"pseudorep"] <- sample(subgroup, size=size, replace=FALSE)
Pseudorep1 <- Maledata.df[Maledata.df$pseudorep %in% 1,]
dim(Pseudorep1)
Pseudorep2 <- Maledata.df[Maledata.df$pseudorep %in% 2,]
dim(Pseudorep2)
Pseudorep3 <- Maledata.df[Maledata.df$pseudorep %in% 3,]
dim(Pseudorep3)

Pseudorep1 <- t(Pseudorep1)
Pseudorep2 <- t(Pseudorep2)
Pseudorep3 <- t(Pseudorep3)

Pseudobulk1 <- rowSums(Pseudorep1)
Pseudobulk2 <- rowSums(Pseudorep2)
Pseudobulk3 <- rowSums(Pseudorep3)
dim(Pseudobulk1)
dim(Pseudobulk2)
dim(Pseudobulk3)
MalePseudoreps <- cbind(Pseudobulk1, Pseudobulk2, Pseudobulk3)
head(MalePseudoreps)



Femaledata <- Female@assays[["RNA"]]@counts
dim(Femaledata)
Femaledata.df <- as.data.frame(t(Femaledata))
dim(Femaledata.df)
subgroupf <- rep(1:3, ceiling(nrow(Femaledata.df)/3))
subgroupf
dim(subgroupf)
Femaledata.df$pseudorep <- 0 
size <- nrow(Femaledata.df)
size
Femaledata.df[,"pseudorep"] <- sample(subgroupf, size=size, replace=FALSE)
Pseudorep1f <- Femaledata.df[Femaledata.df$pseudorep %in% 1,]
Pseudorep2f <- Femaledata.df[Femaledata.df$pseudorep %in% 2,]
Pseudorep3f <- Femaledata.df[Femaledata.df$pseudorep %in% 3,]

Pseudorep1f <- t(Pseudorep1f)
Pseudorep2f <- t(Pseudorep2f)
Pseudorep3f <- t(Pseudorep3f)

Pseudobulk1f <- rowSums(Pseudorep1f)
Pseudobulk2f <- rowSums(Pseudorep2f)
Pseudobulk3f <- rowSums(Pseudorep3f)
dim(Pseudobulk1f)
dim(Pseudobulk2f)
dim(Pseudobulk3f)
FemalePseudoreps <- cbind(Pseudobulk1f, Pseudobulk2f, Pseudobulk3f)
head(FemalePseudoreps)



Combined1 <- cbind(MalePseudoreps, FemalePseudoreps)
Combined2 <- as.matrix(Combined1)
Combined <- round(Combined2)
cond1 <- "Intact"
cond2 <- "Primed"
sampleTable <- data.frame(condition=factor(c(rep(cond1, 3), rep(cond2, 3))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)


#restargeted <- res[filtered.genelist,]
#restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(res, file=paste0(output,i,"_sDEGpadjust.csv"))
res<- na.omit(res)
padjcounts <- nrow(res[res$padj < 0.05,])
df <- data.frame(i, padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- res[res$padj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(output,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(output,"totalsigDEGs.csv"))
genes <- rownames(master_list)
genes <- unlist(genes)
genes <- unique(genes)
write.csv(genes, file=paste0(output, "uniquesigDEGs.csv"))
