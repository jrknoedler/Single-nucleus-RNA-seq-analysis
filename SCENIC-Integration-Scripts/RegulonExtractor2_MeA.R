#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(SCopeLoomR)
library(Seurat)
library(AUCell)

output <- "RegulonsCompiled/MeA/MeA_Unfiltered"

loom <- open_loom("pySCENIC_Singularity/MeA_DefaultNES_Final/auc_mtx_filtered.loom") 
regulons_incidMat <- get_regulons(loom, attrName="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)



MvUDEGS <- read.table("Genelists/MeA_MvU_1.5cutoff.txt")
MvUDEGS <- unlist(MvUDEGS)
MvUDEGS <- as.matrix(MvUDEGS)

PvUDEGS <- read.table("Genelists/MeA_PvU_1.5cutoff.txt")
PvUDEGS <- unlist(PvUDEGS)
PvUDEGS <- as.matrix(PvUDEGS)

MvPDEGS <- read.table("Genelists/MeA_MvP_1.5cutoff.txt")
MvPDEGS <- unlist(MvPDEGS)
MvPDEGS <- as.matrix(MvPDEGS)

AllDegs <- rbind(MvPDEGS,MvUDEGS,PvUDEGS)
AllDegs <- unique(AllDegs)


#sox5 <- data.frame(sox5)
#sox5
#sox5di <- sox5[MvUDEGS,]
#sox5di
regs <- names(regulons)
regs
regs <- unlist(regs)
regs
regs <- data.frame(regs)
regs
write.csv(regs, file=paste0(output,"RegulonNames.csv"))

master_list1 = data.frame()
total_list1=data.frame()
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
total_list1=rbind(total_list1,DEGsReg)
}

final1 <- cbind(master_list1, regs)
write.csv(final1, file=paste0(output,"regulonoverlapMvU.csv"))

master_list2 = data.frame()
total_list2=data.frame()
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
total_list2=rbind(total_list2,DEGsReg)
}

final2 <- cbind(master_list2, regs)
write.csv(final2, file=paste0(output,"regulonoverlapMvP.csv"))

master_list3 = data.frame()
total_list3=data.frame()
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
total_list3=rbind(total_list3,DEGsReg)
}

sRegDegs <- rbind(total_list1,total_list2)
sRegDegs <- unique(sRegDegs)
write.csv(sRegDegs, file=paste0(output,"sDEGsinRegulons.csv"))
eRegDegs <- unique(total_list3)
write.csv(eRegDegs, file=paste0(output,"eDEGsinRegulons.csv"))
TotalDegs <- rbind(sRegDegs,eRegDegs)
TotalDegs <- unique(TotalDegs)
write.csv(TotalDegs, file=paste0(output,"allDEGsinRegulons.csv"))

final3 <- cbind(master_list3, regs)
write.csv(final3, file=paste0(output,"regulonoverlapPvU.csv"))

master_list4 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, AllDegs)
DEGsReg <- data.frame(DEGsReg)
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list4=rbind(master_list4, df2)
}

final4 <- cbind(master_list4, regs)
write.csv(final4, file=paste0(output,"regulonoverlapallVMHDEGs.csv"))

loom <- close_loom()

