#!/usr/bin/env Rscript

output <- "Seurat/VMH_IndependentAnalysis/VMH_Regulons"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(loomR)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
pyScenicLoomFile <- file.path("pySCENIC_Singularity/VMH_DefaultNES_Final/auc_mtx_filtered.loom")
loom <- open_loom(pyScenicLoomFile,mode="r")
exprMat <- get_dgem(loom)
cellInfo <- get_cellAnnotation(loom)
head(cellInfo)
regulons_incidMat <- get_regulons(loom, attrName='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulonsAuc(loom, attrName='RegulonsAUC')
regulonsAUC
regulonsAucThresholds <- get_regulonThresholds(loom)
regulonsAucThresholds
embeddings <- get_embeddings(loom)
embeddings
close_loom(loom)

write.csv(regulonsAucThresholds, file=paste0(output,"AUCthresholds.csv"))

regulonsauc_text <- read.csv("pySCENIC_Singularity/VMH_DefaultNES_Final/auc_mtx_filtered.csv", header=TRUE, check.names=FALSE)
dim(regulonsauc_text)
regulonsauc_text2 <- regulonsauc_text[,-1]
rownames(regulonsauc_text2) <- regulonsauc_text[,1]
regulonsAUC_t <- t(as.matrix(regulonsauc_text2))
regulons_binaryAUC <- ifelse(regulonsAUC_t > regulonsAucThresholds,1,0)
write.csv(regulons_binaryAUC, file=paste0(output,"binarizedmatrix.csv"))
#pdf(file=paste0(output,"binarizetest.pdf"))
#pheatmap(regulons_binaryAUC)
#dev.off()
idents <- mySeurat@active.ident
idents <- data.frame(idents)
regulons_binaryAUC <- t(regulons_binaryAUC)
dim(regulons_binaryAUC)
dim(idents)
head(idents)
merged <- cbind(regulons_binaryAUC, idents)
avg <- merged %>% group_by(idents) %>% summarize_all(.funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames(var="idents")
dim(avg)
avg <- data.frame(avg, check.names=FALSE)
idents <- avg[1]
exp <- avg[c(2:339)]
row.names(exp) <- c(0:27)
exp <- data.frame(exp, check.names=FALSE)
exp <- t(exp)
pdf(file=paste0(output,"binarizedavg.pdf"))
pheatmap(exp,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()



regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
regulonActivity_byCellType_Scaled
pdf(file=paste0(output,"allregulonheatmap.pdf"))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()


Degulators <- read.table("RegulonsCompiled/VMH/Degulators.txt")
Degulators <- unlist(Degulators)
Degulators <- as.matrix(Degulators)
DegulatorsBinarized <- exp[Degulators,]
Degulators <- regulonsAUC[Degulators]



pdf(file=paste0(output,"DegulatorsBinarizedavg.pdf"))
pheatmap(DegulatorsBinarized,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

Degulons <- read.table("RegulonsCompiled/VMH/VMH_SCENICDegulons.txt")
Degulons <- unlist(Degulons)
DegulonsBinarized <- exp[Degulons,]
Degulons <- as.matrix(Degulons)
Degulons <- regulonsAUC[Degulons]


pdf(file=paste0(output,"DegulonsBinarized.pdf"))
pheatmap(DegulonsBinarized,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
DegulatorActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Degulators)[,cells]))
DegulatorActivity_byCellType_Scaled <- t(scale(t(DegulatorActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"10geneDegulatorheatmap.pdf"))
pheatmap::pheatmap(DegulatorActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()




DegulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Degulons)[,cells]))
DegulonActivity_byCellType_Scaled <- t(scale(t(DegulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Degulonheatmap.pdf"))
pheatmap::pheatmap(DegulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

Negulons <- read.table("RegulonsCompiled/VMH/VMH_Negulons.txt")
Negulons <- unlist(Negulons)
Negulons <- as.matrix(Negulons)
NegulonsBinarized <- exp[Negulons,]
Negulons <- regulonsAUC[Negulons]


pdf(file=paste0(output,"NegulonsBinarized.pdf"))
pheatmap(NegulonsBinarized,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

NegulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Negulons)[,cells]))
NegulonActivity_byCellType_Scaled <- t(scale(t(NegulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Negulonheatmap.pdf"))
pheatmap::pheatmap(NegulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

Posulons <- read.table("RegulonsCompiled/VMH/VMH_Posulons.txt")
Posulons <- unlist(Posulons)
Posulons <- as.matrix(Posulons)
PosulonsBinarized <- exp[Posulons,]
Posulons <- regulonsAUC[Posulons]


pdf(file=paste0(output,"PosulonsBinarized.pdf"))
pheatmap(PosulonsBinarized,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()


PosulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Posulons)[,cells]))
PosulonActivity_byCellType_Scaled <- t(scale(t(PosulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Posulonheatmap.pdf"))
pheatmap::pheatmap(PosulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()