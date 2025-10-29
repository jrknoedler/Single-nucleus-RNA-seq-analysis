#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_Regulons"
regoutput <- "RegulonsCompiled/MeA/Sexmarkersperclust/SexmarkerRegulonOverlap"

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
library(stringr)

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_Test_Striatumfiltered2_sexclude_malat1regress_30pcsres1.2FINAL.rds")
pyScenicLoomFile <- file.path("pySCENIC_Singularity/MeA_DefaultNES_Final/auc_mtx_filtered.loom")
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

regs <- names(regulons)
regs <- unlist(regs)
regs <- data.frame(regs)

files <- list.files(path="RegulonsCompiled/MeA/Sexmarkersperclust/", pattern="*.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)
lapply(files, function(x){
markers <- read.csv(x, header=TRUE, row.names=1)
sigmarkers <- markers[markers$p_val_adj < 0.05,]
sig_genes <- rownames(sigmarkers)
sig_genes <- unlist(sig_genes)
sig_genes <- as.matrix(sig_genes)
master_list1 = data.frame()
total_list1=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, sig_genes)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(x,"_regulonoverlap.csv"))
})

files2 <- list.files(path="Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP/", pattern="*sDEGpadjust.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

lapply(files2, function(x){
markers <- read.csv(x, header=TRUE, row.names=1)
markers <- na.omit(markers)
sigmarkers <- markers[markers$padj < 0.05,]
sig_genes <- rownames(sigmarkers)
sig_genes <- unlist(sig_genes)
sig_genes <- as.matrix(sig_genes)
master_list1 = data.frame()
total_list1=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, sig_genes)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(x,"_regulonoverlap.csv"))
})

files3 <- list.files(path="Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvU/", pattern="*sDEGpadjust.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

lapply(files3, function(x){
markers <- read.csv(x, header=TRUE, row.names=1)
markers <- na.omit(markers)
sigmarkers <- markers[markers$padj < 0.05,]
sig_genes <- rownames(sigmarkers)
sig_genes <- unlist(sig_genes)
sig_genes <- as.matrix(sig_genes)
master_list1 = data.frame()
total_list1=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, sig_genes)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(x,"_regulonoverlap.csv"))
})

files4 <- list.files(path="Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/PvU/", pattern="*sDEGpadjust.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

lapply(files4, function(x){
markers <- read.csv(x, header=TRUE, row.names=1)
markers <- na.omit(markers)
sigmarkers <- markers[markers$padj < 0.05,]
sig_genes <- rownames(sigmarkers)
sig_genes <- unlist(sig_genes)
sig_genes <- as.matrix(sig_genes)
master_list1 = data.frame()
total_list1=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, sig_genes)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(x,"_regulonoverlap.csv"))
})

Sexmarkers <- read.table("RegulonsCompiled/MeA/DEGMarkers.txt")
Sexmarkers <- unlist(Sexmarkers)
Sexmarkers <- unique(Sexmarkers)
Sexmarkers <- as.matrix(Sexmarkers)

ClustDEGs <- read.table("RegulonsCompiled/MeA/ClusterSigDEGs.txt")
ClustDEGs <- unlist(ClustDEGs)
ClustDEGs <- unique(ClustDEGs)
ClustDEGs <- as.matrix(ClustDEGs)

master_list1 = data.frame()
total_list1=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, Sexmarkers)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(output,"regulonoverlapSexmarkers.csv"))

master_list2 = data.frame()
total_list2=data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, ClustDEGs)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list2=rbind(master_list2, df2)
total_list2=rbind(total_list2,DEGsReg)
}


final2 <- cbind(master_list2, regs)
write.csv(final2, file=paste0(output,"clustsigDEGSinregulons.csv"))

cellsAUC <- regulonsAUC[c("Nfib(+)", "Tef(+)")]
pdf(file=paste0(output,"Nfibhist.pdf"))
AUCell_plotHist(cellsAUC, aucThr=c(0.23,0.128))
dev.off()

cellsAUC <- regulonsAUC[c("Tcf4(+)")]
pdf(file=paste0(output,"Tcf4hist.pdf"))
AUCell_plotHist(cellsAUC, aucThr=0.246649891)
dev.off()

cellsAUC <- regulonsAUC[c("Tef(+)")]
pdf(file=paste0(output,"Tefhist.pdf"))
AUCell_plotHist(cellsAUC, aucThr=0.128)
dev.off()


write.csv(regulonsAucThresholds, file=paste0(output,"AUCthresholds.csv"))

regulonsauc_text <- read.csv("pySCENIC_Singularity/MeA_DefaultNES_Final/auc_mtx_filtered.csv", header=TRUE, check.names=FALSE)
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
exp <- avg[c(2:201)]
row.names(exp) <- c(0:34)
exp <- data.frame(exp, check.names=FALSE)
exp <- t(exp)
pdf(file=paste0(output,"binarizedavg.pdf"))
pheatmap(exp,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

cols <- colorRampPalette(c("white","purple"))(100)
expcleaned <- ifelse(exp <0.3, 0, exp)
pdf(file=paste0(output,"binarizedavgcleaned.pdf"))
pheatmap(expcleaned, show_rownames=T, border_color=NA, na_col="white", color=magma(2000))
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


Degulators <- read.table("RegulonsCompiled/MeA/AllDegulators.txt")
Degulators <- unlist(Degulators)
Degulators <- as.matrix(Degulators)
DegulatorsBinarized <- exp[Degulators,]
Degulators <- regulonsAUC[Degulators]



pdf(file=paste0(output,"DegulatorsBinarizedavg.pdf"))
pheatmap(DegulatorsBinarized,show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

Degulons <- read.table("RegulonsCompiled/MeA/MeA_SCENICDegulons.txt")
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
                   treeheight_row=10, treeheight_col=10, border_color=NA, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()




DegulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Degulons)[,cells]))
DegulonActivity_byCellType_Scaled <- t(scale(t(DegulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Degulonheatmap.pdf"))
pheatmap::pheatmap(DegulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

Negulons <- read.table("RegulonsCompiled/MeA/MeA_Negulons.txt")
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

Posulons <- read.table("RegulonsCompiled/MeA/MeA_Posulons.txt")
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


#pdf(file=paste0(output,"_allregbindotplot.pdf"), width=16)
#ggplot(exp) + 
#  geom_point() + 
#  scale_color_viridis_c() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  ylab('') +
#  theme(axis.ticks = element_blank())
#dev.off()


df <- data.frame(regulonActivity_byCellType_Scaled)
df <- t(df)
df <- as.matrix(df)
topregs <- colnames(df)[max.col(df, ties.method="first")]
pdf(file=paste0(output,"topRegbyclust.pdf"))
pheatmap(topregs, color=viridis(2000), treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

top5 <- df %>% top_n(n=5)
top5
pdf(file=paste0(output,"top5Regbyclust.pdf"))
pheatmap(top5, color=viridis(2000), treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

df2 <- data.frame(DegulatorActivity_byCellType_Scaled)
top5dr <- df2 %>% top_n(n=5)
top5dr
pdf(file=paste0(output,"top5Degulatorbyclust.pdf"))
pheatmap(top5dr, color=viridis(2000), treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

topregs <- colnames(flipped)[max.col(flipped, ties.method="first")]
topregs
pdf(file=paste0(output,"topRegbyclust.pdf"))
pheatmap(topregs, color=viridis(2000), treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

top5 <- flipped %>% top_n(n=5)
pdf(file=paste0(output,"top5Regbyclust.pdf"))
pheatmap(top5, color=viridis(2000), treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
#rss <- calcRSS(AUC=getAUC(regulonsAUC), cellAnnotation=cellInfo[colnames(regulonsAUC), "seurat_clusters"])
#pdf(file=paste0(output,"RSS.plot"))
#rssPlot <- plotRSS(rss)
#plotly::ggploty(rssPlot$plot)
#dev.off()