#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(SCopeLoomR)
library(Seurat)
library(AUCell)

output <- "RegulonsCompiled/MeA/MeA_"

loom <- open_loom("pySCENIC_Singularity/MeAMerged/auc_mtx_filtered.loom") 
regulons_incidMat <- get_regulons(loom, attrName="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
names(regulons)
regulons$AR
regulons$Ar
names(regulons$Foxo1)

arg <- regulons$Egr2
write.csv(arg, file=paste0(output,"arg1.csv"))
arg <- regulons$Egr3
write.csv(arg, file=paste0(output,"arg2.csv"))
arg <- regulons$Fosb
write.csv(arg, file=paste0(output,"arg3.csv"))
arg <- regulons$Jun
write.csv(arg, file=paste0(output,"arg4.csv"))
arg <- regulons$Junb
write.csv(arg, file=paste0(output,"arg5.csv"))
arg <- regulons$Nfia
write.csv(arg, file=paste0(output,"arg6.csv"))
arg <- regulons$Nfix
write.csv(arg, file=paste0(output,"arg7.csv"))
arg <- regulons$Npas2
write.csv(arg, file=paste0(output,"arg8.csv"))
arg <- regulons$Nr3c2
write.csv(arg, file=paste0(output,"arg9.csv"))
arg <- regulons$Nr4a3
write.csv(arg, file=paste0(output,"arg10.csv"))
arg <- regulons$Rorb
write.csv(arg, file=paste0(output,"arg11.csv"))
arg <- regulons$Sox5
write.csv(arg, file=paste0(output,"arg12.csv"))

MvUDEGS <- read.table("Genelists/MeA_MvU_1.5cutoff.txt")
MvUDEGS <- unlist(MvUDEGS)
MvUDEGS <- as.matrix(MvUDEGS)

PvUDEGS <- read.table("Genelists/MeA_PvU_1.5cutoff.txt")
PvUDEGS <- unlist(PvUDEGS)
PvUDEGS <- as.matrix(PvUDEGS)

MvPDEGS <- read.table("Genelists/MeA_MvP_1.5cutoff.txt")
PvPDEGS <- unlist(PvUDEGS)
PvUDEGS <- as.matrix(PvUDEGS)



#sox5 <- data.frame(sox5)
#sox5
#sox5di <- sox5[MvUDEGS,]
#sox5di
regs <- names(regulons)
regs <- unlist(regs)
regs <- data.frame(regs)
write.csv(regs, file=paste0(output,"RegulonNames.csv"))

master_list1 = data.frame()
total_list=data.frame()
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
write.csv(final1, file=paste0(output,"POAregulonoverlapMvU.csv"))

master_list2 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, MvPDEGS)
DEGsReg <- data.frame(DEGsReg)
write.csv(DEGsReg, file=paste0(output,i,"MvPDEGs.csv"))
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list2, df2)
}

final2 <- cbind(master_list2, regs)
write.csv(final2, file=paste0(output,"POAregulonoverlapMvP.csv"))

master_list3 = data.frame()
for (i in regulons){
df <- data.frame(i)
total <- nrow(df)
df <- as.matrix(df)
regname <- colnames(i)
DEGsReg <- intersect(df, PvUDEGS)
DEGsReg <- data.frame(DEGsReg)
write.csv(DEGsReg, file=paste0(output,i,"PvUDEGs.csv"))
frac <- nrow(DEGsReg)
df2 <- data.frame(frac, total)
master_list1=rbind(master_list3, df2)
}

final3 <- cbind(master_list3, regs)
write.csv(final3, file=paste0(output,"regulonoverlapPvU.csv"))



loom <- close_loom()

