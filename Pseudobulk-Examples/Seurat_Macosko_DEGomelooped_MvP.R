#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Macosko_KnoedlerDEGs/MvP/MvP_Macosko1.5fold_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(DESeq2)
library(Matrix)
library(matrixStats)
library(tibble)


mySeurat <- readRDS("Seurat/BNSTMacosko/Macosko_BNST_Exploratory_annotated.rds")


idents <- mySeurat@active.ident



genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
total_dimorphic = data.frame()
master_list = data.frame()
DefaultAssay(mySeurat) <- "RNA"

total_dimorphic = data.frame()
master_list = data.frame()
replicates <- c("MALE1","MALE2","MALE3","MALE4","MALE5","MALE6","MALE7","MALE8","FEMALE1","FEMALE2","FEMALE3","FEMALE4","FEMALE5","FEMALE6","FEMALE7")

for (i in 0:50){
try({
ident <- i

cluster <- subset(x=mySeurat, idents=c(i))
count.table = data.frame(dummy=1)
for (g in replicates){

col <- subset(cluster, subset=Replicate==g)
coldata <- col@assays[["RNA"]]@counts
coldata.df <- as.data.frame(coldata)
Pseudobulk <- rowSums(coldata.df)
count.table=cbind(count.table, Pseudobulk)
}
count.table <- count.table[,2:ncol(count.table)]
count.table <- t(count.table)

cond1 <- "Male"
cond2 <- "Female"
sampleTable <- data.frame(condition=factor(c(rep(cond1, 8), rep(cond2, 7))))
rownames(sampleTable) <- colnames(count.table)
dds=DESeqDataSetFromMatrix(count.table, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
restargeted <- res[filtered.genelist,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(restargeted, file=paste0(output,i,"_sDEGpadjust.csv"))
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
