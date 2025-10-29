output <- "Seurat/VMH_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/Top1_MvP"

library(tidyverse)
library(pheatmap)
library(viridis)
library(RColorBrewer)

files <- list.files(path="Seurat/VMH_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP_Sig/", pattern="*.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

master_list = data.frame()
for (i in files) 
try({
DEGs <- read.csv(i, header=TRUE, row.names=1)
DEGs <- na.omit(DEGs)
Sig <- DEGs[DEGs$padj < 0.05,]
Sig <- data.frame(Sig)
Top5 <- Sig %>% top_n(1, abs(log2FoldChange))
genes <- rownames(Top5)
genes <- unlist(genes)
genes <- data.frame(genes)
master_list=rbind(master_list,genes)
})
master_list <- unlist(master_list)
master_list <- unique(master_list)
write.table(master_list, file=paste0(output,"MvP.txt"))
warnings()
df <- data.frame(master_list)
df %>% remove_rownames %>% column_to_rownames(var="master_list")
ord.df<-cbind(rownames(df)[order(rownames(df))], df[order(rownames(df)),])
ord.df
#for (i in files)

#try({
DEGs <- read.csv("Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP_Sig/MvP_DEGome1.5fold_reclustered0_sDEGpadjust.csv", header=TRUE, row.names=1)
DEGs
DEGs <- data.frame(DEGs)
DEGs <- na.omit(DEGs)
Genes <- DEGs[master_list,]
Genes
Genes <- data.frame(Genes)
ordered.Genes<-cbind(rownames(Genes)[order(rownames(Genes))], Genes[order(rownames(Genes)),])
ordered.Genes
#list <- rownames(Genes)
#list
#list <- unlist(list)
#list
#list <- data.frame(list)
Genes <- Genes[,c(3)]
ord.df
ord.df=cbind(ord.df, Genes)
#})
ord.df
write.csv(df, file=paste0(output,"fcgenes.csv"))

df[is.na(data)]=0
paletteLength <- 100
myColor <- colorRampPalette(c("blue","white","red"))(paletteLength)
myBreaks <- c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))
pdf(file=paste0(output,"BNST_MvP.pdf", height=20))
pheatmap(data, color=myColor, breaks=myBreaks, border_color=NA, show_rownames=TRUE)
dev.off()



