library(EnhancedVolcano)
data1 <- read.table("BNST_MvP_Volcano.txt", header=TRUE) 
data2 <- read.table("BNST_MvU_Volcano.txt", header=TRUE)
data3 - read.table("BNST_PvU_Volcano.txt", header=TRUE)

keyvals1 <- ifelse(data1$log2FoldChange < -0.584962501 && data1$padj < 0.05, 'red', ifelse(data1$log2FoldChange > 0.584962501 && data1$padj < 0.05, 'blue', 'gray')) 
keyvals1[is.na(keyvals1)] <- 'gray'
names(keyvals1)[keyvals1=='red'] <- 'Fr' 
names(keyvals1)[keyvals1=='blue'] <- 'M' 
names(keyvals1)[keyvals1=='gray'] <- 'NC' 

keyvals2 <- ifelse(data2$log2FoldChange < -0.584962501 && data2$padj < 0.05, 'green', ifelse(data2$log2FoldChange > 0.584962501 && data2$padj < 0.05, 'blue', 'gray')) 
keyvals2[is.na(keyvals2)] <- 'gray'
names(keyvals2)[keyvals2=='green'] <- 'Fu' 
names(keyvals2)[keyvals2=='blue'] <- 'M' 
names(keyvals2)[keyvals2=='gray'] <- 'NC' 


keyvals3 <- ifelse(data3$log2FoldChange < -0.584962501 && data3$padj < 0.05, 'green', ifelse(data2$log2FoldChange > 0.584962501 && data2$padj < 0.05, 'red', 'gray')) 
keyvals3[is.na(keyvals3)] <- 'gray'
names(keyvals3)[keyvals3=='green'] <- 'Fu' 
names(keyvals3)[keyvals3=='blue'] <- 'M' 
names(keyvals3)[keyvals3=='gray'] <- 'NC' 

pdf(file="BNST_MvP_volcano.pdf")
EnhancedVolcano(data1, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals2, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()

pdf(file="BNST_MvU_volcano.pdf")
EnhancedVolcano(data1, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals2, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()

pdf(file="BNST_PvU_volcano.pdf")
EnhancedVolcano(data3, lab = NA, x='log2FoldChange', y='padj', gridlines.major=FALSE, gridlines.minor=FALSE, colCustom=keyvals3, selectLab=c(), pCutoff=NA, FCcutoff=0.584962501)
dev.off()