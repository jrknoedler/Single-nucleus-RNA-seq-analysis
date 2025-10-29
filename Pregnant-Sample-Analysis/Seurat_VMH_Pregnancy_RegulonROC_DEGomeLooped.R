#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
pregvpoutput <- "Seurat/VMH_Barcode_Transfer/Cluster_DEGs/Regulon_ROC/Pregnant_vs_Primed/PregvPrimed_Regulons_"
pregvuoutput <- "Seurat/VMH_Barcode_Transfer/Cluster_DEGs/Regulon_ROC/Pregnant_vs_Unprimed/PregvUnprimed_Regulons_"

library(Seurat)
library(reticulate)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/Pregnancy_Full_Dataset_Analysis/Merged_Pregnancy_with_Positive_RegulonsVMH_withPosRegulons.rds")
DefaultAssay(mySeurat) <- "RNA"

clustno <- nlevels(Idents(mySeurat))
new.clust.ids <- c(1:clustno)

names(new.clust.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.clust.ids)

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone

BNST <- read.csv("Seurat/Pregnancy_Full_Dataset_Analysis/Regulons_1stpass/Regulons_OnlyPos_lowNESregulonoverlap_BNSTpregDEGs.csv", header=TRUE, row.names=1)

#Get relevant DEGs
PregvEst <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_PvFr_All1.5.txt")
PregvEst <- unlist(PregvEst)
PregvDiest <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_PregvFnr_All1.5.txt")
PregvDiest <- unlist(PregvDiest)
PubDEGs <- read.table("topGO/Total_ByRegion/VMH_Genesonly.txt")
PubDEGs <- unlist(PubDEGs)
DEG_int <- union(PregvEst, PregvDiest)
DEGs_Master <- union(DEG_int, PubDEGs)
VMHPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_TotalUpPreg.txt")
VMHPregup <- unlist(VMHPregup)


#DefaultAssay(VMHcells) <- "Regulons"
#VMH_reg <- FindAllMarkers(VMHcells, test.use="roc", verbose=TRUE, logfc.threshold=0, return.thresh=0)

##Get all regulon names and associated genes
AllRegulons <- BNST[,1]
AllRegulons_ed <- gsub("\\(\\+\\)","",AllRegulons)

##Subset regulons that match genes expressed in >25% of cells in at least one cluster, then cry when you think about wasted analysis time
Pct <- DotPlot(mySeurat, features=AllRegulons_ed)
data <- Pct$data
head(data)
data.df <- data.frame(data)
exp.regs.filt <- data.df[data.df$pct.exp >= 40,]
exp.regs <- exp.regs.filt$features.plot
exp.regs <- unique(exp.regs)
Exp.DEG.regs <- intersect(exp.regs, DEGs_Master)
Exp.Pregup.DEGs <- intersect(exp.regs, VMHPregup)
exp.regs.Regulons <- paste0(exp.regs,"(+)")

mySeurat$celltype.Hormone <- paste(Idents(mySeurat), mySeurat$Hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.Hormone"
total_dimorphic = data.frame()
master_list = data.frame()
DefaultAssay(mySeurat) <- "Regulons"
for (i in 1:clustno){
try({
ident1 <- paste0(i,"_Pregnant")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(mySeurat, ident.1=ident1, ident.2=ident2, test.use="roc", verbose=TRUE, logfc.threshold=0, return.thresh=0, features=exp.regs.Regulons)
write.csv(sex.dimorphism, file=paste0(pregvpoutput,i,"_Regulons.csv"))
})
}

for (i in 1:clustno){
try({
ident1 <- paste0(i,"_Pregnant")
ident2 <- paste0(i,"_Unprimed")
sex.dimorphism <- FindMarkers(mySeurat, ident.1=ident1, ident.2=ident2, test.use="roc", verbose=TRUE, logfc.threshold=0, return.thresh=0, features=exp.regs.Regulons)
write.csv(sex.dimorphism, file=paste0(pregvuoutput,i,"_Regulons.csv"))

})
}

Idents(mySeurat) <- "Hormone"
PregvPrimed <- FindMarkers(mySeurat, ident.1="Pregnant", ident.2="Primed", test.use="roc", verbose=TRUE, logfc.threshold=0, return.thresh=0)
write.csv(PregvPrimed, file=paste0(pregvpoutput,"PregvPrimed_TotalRegs.csv"))
PregvUnprimed <- FindMarkers(mySeurat, ident.1="Pregnant", ident.2="Unprimed", test.use="roc", verbose=TRUE, logfc.threshold=0, return.thresh=0)
write.csv(PregvUnprimed, file=paste0(pregvuoutput,"PregvUnprimed_TotalRegs.csv"))