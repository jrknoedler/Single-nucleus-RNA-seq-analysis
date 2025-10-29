library(dplyr)
library(tidyverse)
output <- "DEclustcounts/BNST/"



files <- list.files(path="Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/PvU_Sig/", pattern="*padjust.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)


master_list = data.frame()
for (i in files) 

try({
DEGs <- read.csv(i, header=TRUE, row.names=1)
DEGs.df <- as.data.frame(DEGs)
DEGs.df
DEGs.df <- na.omit(DEGs.df)
DEGs.df
DEGs.df <- DEGs.df[DEGs.df$padj <0.05,]
DEGs.df
DEGs.df <- rownames_to_column(DEGs.df, var="gene")
DEGs.df
siggenes <- DEGs.df["gene"]
master_list <- rbind(master_list, siggenes)


})
master_list
total <- table(master_list$gene)
total
write.csv(total, file=paste0(output,"PvU_clustsigcounts.csv"))