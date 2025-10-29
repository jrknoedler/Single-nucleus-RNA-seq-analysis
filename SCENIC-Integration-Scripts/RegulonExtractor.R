#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(SCopeLoomR)
library(Seurat)
library(AUCell)

output <- "RegulonsCompiled/VMH_DefaultNESfinal_"

loom <- open_loom("pySCENIC_Singularity/VMH_DefaultNES_Final/auc_mtx_filtered.loom") 
regulons_incidMat <- get_regulons(loom, attrName="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
names(regulons)
regulons$AR
regulons$Ar
regulons$Foxo1
MvUDEGS <- read.table("Genelists/VMH_MvU_1.5.txt")
MvUDEGS <- unlist(MvUDEGS)
MvUDEGS <- as.matrix(MvUDEGS)

PvUDEGS <- read.table("Genelists/VMH_PvU_1.5.txt")
PvUDEGS <- unlist(PvUDEGS)
PvUDEGS <- as.matrix(PvUDEGS)

MvPDEGS <- read.table("Genelists/VMH_MvP_1.5.txt")
MvPDEGS <- unlist(MvPDEGS)
MvUDEGS <- as.matrix(MvPDEGS)



#sox5 <- data.frame(sox5)
#sox5
#sox5di <- sox5[MvUDEGS,]
#sox5di
regs <- names(regulons)
regs <- unlist(regs)
regs <- data.frame(regs)
write.csv(regs, file=paste0(output,"RegulonNames.csv"))

master_list1 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, MvUDEGS)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list1, df2)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(output,"regulonoverlapMvU.csv"))

master_list2 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, MvPDEGS)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list2=rbind(master_list2, df2)
}

final2 <- cbind(master_list2, regs)
write.csv(final2, file=paste0(output,"regulonoverlapMvP.csv"))

master_list3 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, PvUDEGS)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list3=rbind(master_list3, df2)
}

final3 <- cbind(master_list3, regs)
write.csv(final3, file=paste0(output,"regulonoverlapPvU.csv"))



loom <- close_loom()

