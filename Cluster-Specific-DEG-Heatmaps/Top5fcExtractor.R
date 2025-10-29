output <- "Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/Top5_"

library(tidyverse)

files <- list.files(path="Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/PvU_Sig/", pattern="*.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)

master_list = data.frame()
for (i in files) 
try({
DEGs <- read.csv(i, header=TRUE, row.names=1)
DEGs <- na.omit(DEGs)
Sig <- DEGs[DEGs$padj < 0.05,]
Sig <- data.frame(Sig)
Top5 <- Sig %>% top_n(5, abs(log2FoldChange))
genes <- rownames(Top5)
genes <- unlist(genes)
genes <- data.frame(genes)
master_list=rbind(master_list,genes)
})

master_list <- unique(master_list)
write.csv(master_list, file=paste0(output,"PvU.txt"))
warnings()



