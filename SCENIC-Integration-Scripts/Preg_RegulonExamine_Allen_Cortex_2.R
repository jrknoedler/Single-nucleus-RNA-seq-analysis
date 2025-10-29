#!/usr/bin/env Rscript

output <- "Seurat/VMH_Barcode_Transfer/CortexRegulons_1stpass"

library(Seurat)
#library(pheatmap)
#library(viridis)
#library(RColorBrewer)
#library(dplyr)
#library(tidyverse)
#library(SCENIC)
library(Matrix)
#library(AUCell)
#library(loomR)
library(future)


options(future.globals.maxSize=8000*1024^2)
path="/scratch/users/knoedler/VMH_Allen/Tests/"

VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
DefaultAssay(VMH) <- "RNA"
Allen.data <- read.csv("Allen_Data/Subsampled_30K.csv", header=TRUE, row.names=1)
Allen.data <- t(Allen.data)
Allen <- CreateSeuratObject(counts=Allen.data, project="Allen_10x", min.cells=10, min.features=200)

VMH$Region <- "VMH"
Allen$Region <- "Cortex"

mySeurat <- merge(x=VMH, y=Allen, add.cell.ids=c("VMH","Allen"))

regulons.master = data.frame()
dim(regulons.master)
dirs=list.dirs(path, recursive=FALSE)
dirs
for(i in dirs){
try({
regulons <- read.csv(paste0(i,"/auc_mtx_filtered.csv"), check.names=FALSE, header=TRUE, row.names=1)
regulons <- t(regulons)
regulons.master = rbind(regulons.master, regulons)
})
}
dim(regulons.master)
regs <- as.matrix(regulons.master)
mySeurat[["Regulons"]] <- CreateAssayObject(data = regs)
mySeurat

VMH2 <- subset(mySeurat, subset=Region=="VMH")
saveRDS(VMH2, file=paste0(output,"VMH.rds"))

