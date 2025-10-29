output <- "Seurat/Heatmaps/POA/MvP/"



files <- list.files(path="Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MvP_Sig/", pattern="*.csv", full.names=TRUE, include.dirs=FALSE, recursive=FALSE)
sigdegs <- read.table("Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/AllMvPsig.txt", header=FALSE)
sigdegs <- unlist(sigdegs)
sigdegs <- as.matrix(sigdegs)

master_list = data.frame()
for (i in files) 

try({
DEGs <- read.csv(i, header=TRUE, row.names=1)
DEGs.df <- as.data.frame(DEGs)
Sig <- DEGs.df[sigdegs,]
Sig
Sig <- na.omit(Sig)
Sig <- Sig[c("log2FoldChange","padj")]
Sig
Sig2 <- replace(Sig$log2FoldChange, Sig$padj >=0.05, NA)
Sig2
Sig3 <- cbind(Sig, Sig2)
Sig3
write.csv(Sig3, file=paste0(i, "sigfc.csv"))
})