#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(SCopeLoomR)
library(Seurat)
library(AUCell)

output <- "RegulonsCompiled/BNST/BNST_"

loom <- open_loom("pySCENIC_Singularity/BNSTMerged_lowNES/auc_mtx_filtered.loom") 
regulons_incidMat <- get_regulons(loom, attrName="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
names(regulons)
Arx <- data.frame(regulons$Arx)
write.csv(Arx, file=paste0(output,"arx.csv"))
regulons$Ar
names(regulons$Foxo1)
Degulons <- read.table("RegulonsCompiled/BNST/BNSTDegulons.txt")

Degulons <- unlist(Degulons)


arg <- regulons$Arx
write.csv(arg, file=paste0(output,"arg1.csv"))
arg <- regulons$Atf2
write.csv(arg, file=paste0(output,"arg2.csv"))
arg <- regulons$Bbx
write.csv(arg, file=paste0(output,"arg3.csv"))
arg <- regulons$Bcl6
write.csv(arg, file=paste0(output,"arg4.csv"))
arg <- regulons$Bdp1
write.csv(arg, file=paste0(output,"arg5.csv"))
arg <- regulons$Ccnt2
write.csv(arg, file=paste0(output,"arg6.csv"))
arg <- regulons$Chd1
write.csv(arg, file=paste0(output,"arg7.csv"))
arg <- regulons$Clock
write.csv(arg, file=paste0(output,"arg8.csv"))
arg <- regulons$Creb3l1
write.csv(arg, file=paste0(output,"arg9.csv"))
arg <- regulons$Egr1
write.csv(arg, file=paste0(output,"arg1.csv"))
arg <- regulons$Egr3
write.csv(arg, file=paste0(output,"arg10.csv"))
arg <- regulons$Ep300
write.csv(arg, file=paste0(output,"arg11.csv"))
arg <- regulons$Esrrg
write.csv(arg, file=paste0(output,"arg12.csv"))
arg <- regulons$Fos
write.csv(arg, file=paste0(output,"arg13.csv"))
arg <- regulons$Fosb
write.csv(arg, file=paste0(output,"arg14.csv"))
arg <- regulons$Fosl2
write.csv(arg, file=paste0(output,"arg15.csv"))
arg <- regulons$Foxj2
write.csv(arg, file=paste0(output,"arg16.csv"))
arg <- regulons$Foxp2
write.csv(arg, file=paste0(output,"arg17.csv"))
arg <- regulons$Glis2
write.csv(arg, file=paste0(output,"arg18.csv"))
arg <- regulons$Gtf2ird1
write.csv(arg, file=paste0(output,"arg19.csv"))
arg <- regulons$Hivep1
write.csv(arg, file=paste0(output,"arg20.csv"))
arg <- regulons$Hsf1
write.csv(arg, file=paste0(output,"arg21.csv"))
arg <- regulons$Jun
write.csv(arg, file=paste0(output,"arg22.csv"))
arg <- regulons$Kdm5a
write.csv(arg, file=paste0(output,"arg23.csv"))
arg <- regulons$Kdm5b
write.csv(arg, file=paste0(output,"arg24.csv"))
arg <- regulons$Mef2a
write.csv(arg, file=paste0(output,"arg25.csv"))
arg <- regulons$Meis3
write.csv(arg, file=paste0(output,"arg26.csv"))
arg <- regulons$Mga
write.csv(arg, file=paste0(output,"arg27.csv"))
arg <- regulons$Nfia
write.csv(arg, file=paste0(output,"arg28.csv"))
arg <- regulons$Nfix
write.csv(arg, file=paste0(output,"arg29.csv"))
arg <- regulons$Npas2
write.csv(arg, file=paste0(output,"arg30.csv"))
arg <- regulons$Otx2
write.csv(arg, file=paste0(output,"arg31.csv"))
arg <- regulons$Pgr
write.csv(arg, file=paste0(output,"arg32.csv"))
arg <- regulons$Plagl1
write.csv(arg, file=paste0(output,"arg33.csv"))
arg <- regulons$Ppargc1a
write.csv(arg, file=paste0(output,"arg34.csv"))
arg <- regulons$Prox1
write.csv(arg, file=paste0(output,"arg35.csv"))
arg <- regulons$Pura
write.csv(arg, file=paste0(output,"arg36.csv"))
arg <- regulons$Rfx3
write.csv(arg, file=paste0(output,"arg37.csv"))
arg <- regulons$Sall3
write.csv(arg, file=paste0(output,"arg38.csv"))
arg <- regulons$Sox5
write.csv(arg, file=paste0(output,"arg39.csv"))
arg <- regulons$Sox6
write.csv(arg, file=paste0(output,"arg40.csv"))
arg <- regulons$Tcf4
write.csv(arg, file=paste0(output,"arg41.csv"))
arg <- regulons$Tfdp2
write.csv(arg, file=paste0(output,"arg42.csv"))
arg <- regulons$Thrb
write.csv(arg, file=paste0(output,"arg43.csv"))
arg <- regulons$Trps1
write.csv(arg, file=paste0(output,"arg44.csv"))
arg <- regulons$Zfp445
write.csv(arg, file=paste0(output,"arg45.csv"))
arg <- regulons$Zfp667
write.csv(arg, file=paste0(output,"arg46.csv"))
arg <- regulons$Zic1
write.csv(arg, file=paste0(output,"arg47.csv"))




degulonlist = data.frame()
for (g in Degulons) {
degulon <- data.frame(regulons$g)
#degulon <- data.frame(degulon)
degulonlist =rbind(degulonlist, degulon)
}
write.csv(degulonlist, file=paste0(output,"degulontargets.csv"))


MvUDEGS <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
MvUDEGS <- unlist(MvUDEGS)
MvUDEGS <- as.matrix(MvUDEGS)

PvUDEGS <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
PvUDEGS <- unlist(PvUDEGS)
PvUDEGS <- as.matrix(PvUDEGS)

MvPDEGS <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
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

