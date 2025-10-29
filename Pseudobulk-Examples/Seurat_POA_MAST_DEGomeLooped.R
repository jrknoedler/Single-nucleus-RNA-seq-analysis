#!/usr/bin/env Rscript

pregvpoutput <- "Seurat/POA_PregIntegrated/Cluster_DEGs/MAST_noFilter/Preg_v_Primed/PregvPrimed_DEGome_"
pregvuoutput <- "Seurat/POA_PregIntegrated/Cluster_DEGs/MAST_noFilter/Preg_v_Unprimed/PregvUnprimed_DEGome_"

library(Seurat)
library(reticulate)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")

clustno <- nlevels(Idents(mySeurat))
new.clust.ids <- c(1:clustno)

names(new.clust.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.clust.ids)

genelist <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_PregvFrr_All1.5.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.PregvPrimed <- intersect(genelist, genes.10x)
mySeurat$celltype.Hormone <- paste(Idents(mySeurat), mySeurat$Hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.Hormone"
total_dimorphic = data.frame()
master_list = data.frame()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
for (i in 1:clustno){
try({
ident1 <- paste0(i,"_Pregnant")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(mySeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.PregvPrimed, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(pregvpoutput,i,"_MASTDEGOME.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
padjcounts <- nrow(sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,])
df <- data.frame(i,padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(pregvpoutput,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(pregvpoutput,"totalsigDEGs.csv"))
genes <- rownames(master_list)
genes <-  master_list$gene
genes <- unique(genes)
write.csv(genes, file=paste0(pregvpoutput, "uniquesigDEGs.csv"))


genelist <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_PregvFnr_All1.5.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.PregvUnprimed <- intersect(genelist, genes.10x)
total_dimorphic = data.frame()
master_list = data.frame()
for (i in 1:clustno){
try({
ident1 <- paste0(i,"_Pregnant")
ident2 <- paste0(i,"_Unprimed")
sex.dimorphism <- FindMarkers(mySeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE,features=filtered.PregvUnprimed, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(pregvuoutput,i,"_MASTDEGOME.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
padjcounts <- nrow(sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,])
df <- data.frame(i,padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
clust <- sex.dimorphism[sex.dimorphism$p_val_adj < 0.05,]
clust$gene <- row.names(clust)
master_list = rbind(master_list, clust)
})
}
write.csv(total_dimorphic, file=paste0(pregvuoutput,"_TRAPsDEGsdibyclust.csv"))
write.csv(master_list, file=paste0(pregvuoutput,"totalsigDEGs.csv"))
genes <- master_list$gene
genes <- unlist(genes)
genes <- unique(genes)
write.csv(genes, file=paste0(pregvuoutput, "uniquesigDEGs.csv"))


Idents(mySeurat) <- "Hormone"
PregvPrimed <- FindMarkers(mySeurat, assay="RNA", ident.1 = "Pregnant", ident.2="Primed", min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.PregvPrimed, verbose=TRUE)
write.csv(PregvPrimed, file=paste0(pregvpoutput,"PregvPrimed_Metadegs.csv"))
PregvUnprimed <- FindMarkers(mySeurat, assay="RNA", ident.1 = "Pregnant", ident.2="Unprimed", min.cells.group=0, logfc.threshold=0, min.pct=0, only.pos=FALSE, features=filtered.PregvUnprimed, verbose=TRUE)
write.csv(PregvUnprimed, file=paste0(pregvuoutput,"PregvUnprimed_Metadegs.csv"))