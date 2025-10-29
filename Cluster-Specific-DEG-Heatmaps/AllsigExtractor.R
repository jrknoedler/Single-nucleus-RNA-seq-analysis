output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/AllsigMvP_"

library(tidyverse)

files <- list.files(path="Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP_Sig/", pattern="*.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

master_list = data.frame()
for (i in files) 
try({
DEGs <- read.csv(i, header=TRUE, row.names=1)
DEGs <- na.omit(DEGs)
Sig <- DEGs[DEGs$padj < 0.05,]
genes <- rownames(Sig)
genes <- unlist(genes)
genes <- data.frame(genes)
master_list=rbind(master_list,genes)
})

master_list <- unique(master_list)
write.csv(master_list, file=paste0(output,"MvP.txt"))
warnings()



