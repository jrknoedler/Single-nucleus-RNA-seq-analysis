library(tidyverse)

DEGs <- read.csv("Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP_Sig/MvP_DEGome1.5fold_reclustered0_sDEGpadjust.csv", header=TRUE, row.names=1)
DEGs <- na.omit(DEGs)
DEGs
Sig <- DEGs[DEGs$padj < 0.05,]
Sig
Sig <- data.frame(Sig)
Sig
Top5 <- Sig %>% top_n(5, abs(log2FoldChange))
Top5
Genes <- rownames(Top5)
Genes