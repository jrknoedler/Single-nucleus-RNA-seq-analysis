#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(AUCell)
library(SCopeLoomR)
library(Matrix)

#Set output directory
output <- "C:/Users/Joe/OneDrive - Stanford/Pregnancy_papers/New_Regulons/Pos_neg/"

#Get relevant gene lists and save as matrices


#Only pregnancy-upregulated genes
BNSTDEGs <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_AllvFemale1.5_3.txt")
BNSTDEGs <- unlist(BNSTDEGs)
MeADEGs <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_AllvFemale1.5_2.txt")
MeADEGs <- unlist(MeADEGs)
POADEGs <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_AllvFemale1.5.txt")
POADEGs <- unlist(POADEGs)
VMHDEGs <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_AllvFemale1.5.txt")
VMHDEGs <- unlist(VMHDEGs)
DEGs.master1 <- union(BNSTDEGs,MeADEGs)
DEGs.master2 <- union(DEGs.master1, POADEGs)
DEGs.master3 <- union(DEGs.master2, VMHDEGs)
DEGs.master <- as.matrix(DEGs.master3)

#eDEGs, which we should have done before
BNSTeDEGs <- read.table("Pregnancy Paper/Data/Genelists/BNST_PvU_1.5cutoff.txt")
BNSTeDEGs <- unlist(BNSTeDEGs)
BNSTeDEGs <- as.matrix(BNSTeDEGs)
MeAeDEGs <- read.table("Pregnancy Paper/Data/Genelists/MeA_PvU_1.5cutoff.txt")
MeAeDEGs <- unlist(MeAeDEGs)
MeAeDEGs <- as.matrix(MeAeDEGs)
POAeDEGs <- read.table("Pregnancy Paper/Data/Genelists/POA_PvU_1.5cutoff.txt")
POAeDEGs <- unlist(POAeDEGs)
POAeDEGs <- as.matrix(POAeDEGs)
VMHeDEGs <- read.table("Pregnancy Paper/Data/Genelists/VMH_PvU_1.5.txt")
VMHeDEGs <- unlist(VMHeDEGs)
VMHeDEGs <- as.matrix(VMHeDEGs)
eDEGs.master1 <- union(BNSTeDEGs, MeAeDEGs)
eDEGs.master2 <- union(eDEGs.master1, POAeDEGs)
eDEGs.master3 <- union(eDEGs.master2, VMHeDEGs)
eDEGs.master <- as.matrix(eDEGs.master3)


BNSTBackground <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_Background.txt")
BNSTBackground <- unlist(BNSTBackground)
MeABackground <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_Background.txt")
MeABackground <- unlist(MeABackground)
POABackground <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_Background.txt")
POABackground <- unlist(POABackground)
VMHBackground <- read.table("Pregnancy Paper/Data/Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_Background.txt")
VMHBackground <- unlist(VMHBackground)
Background.master1 <- union(BNSTBackground,MeABackground)
Background.master2 <- union(Background.master1,POABackground)
Background.master <- union(Background.master2, VMHBackground)
Background.master <- as.matrix(Background.master)
head(Background.master)
dim(Background.master)
univt <- nrow(Background.master)

BNSTDEGs <- as.matrix(BNSTDEGs)
head(BNSTDEGs)
dim(BNSTDEGs)
BNSTBackground <- as.matrix(BNSTBackground)
univb <- nrow(BNSTBackground)


MeADEGs <- as.matrix(MeADEGs)

MeABackground <- as.matrix(MeABackground)
univm <- nrow(MeABackground)


POADEGs <- as.matrix(POADEGs)

POABackground <- as.matrix(POABackground)
univp <- nrow(POABackground)


VMHDEGs <- as.matrix(VMHDEGs)

VMHBackground <- as.matrix(VMHBackground)
univv <- nrow(VMHBackground)


BNST.master =data.frame()
MeA.master = data.frame()
POA.master = data.frame()
VMH.master = data.frame()
Total.master = data.frame()
BNSTe.master = data.frame()
MeAe.master = data.frame()
POAe.master = data.frame()
VMHe.master = data.frame()
Totale.master=data.frame()


Total.BNST = data.frame()
Total.MeA = data.frame()
Total.POA = data.frame()
Total.VMH = data.frame()
Total.DEGs = data.frame()
Total.eDEGs = data.frame()
Total.BNSTe = data.frame()
Total.MeAe = data.frame()
Total.POAe = data.frame()
Total.VMHe = data.frame()




loom <- open_loom("X:/Joe_Sherlock_Backup3/Joe_Sherlock_Backup/pySCENIC_Singularity/PseudobulkGRN/GRN_Output/auc_mtx_filtered_2NES_Pseudobulkregulons.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

loom <- close_loom()
###total DEGlists
regs <- names(regulons)
for (i in regs){
try({
    BNST <- data.frame(regulons[i])
    BNST <- as.matrix(BNST)
    regname <- i
    BNSTDEGsReg <- intersect(BNST, BNSTDEGs)
    BNSTDEGsReg <- data.frame(BNSTDEGsReg)
    Total.BNST = rbind(Total.BNST, BNSTDEGsReg)
    BNSTdect <- intersect(BNST,BNSTBackground)
    BNSTdect <- data.frame(BNSTdect)
    BNSTfrac <- nrow(BNSTDEGsReg)
    BNSTrep <- nrow(BNSTDEGsReg)/nrow(BNST)
    total <- nrow(BNSTdect)
    unifrac <- total/univb
    BNSTdf2 <- data.frame(BNSTfrac, total)
    BNSTdf2$pct <- BNSTdf2$BNSTfrac/BNSTdf2$total
    BNSTdf2$pct.degs <- BNSTfrac/nrow(BNSTDEGs)
    BNSTdf3 <- cbind(BNSTdf2, unifrac)
    q1 <- BNSTfrac
    q2 <- (total - BNSTfrac)
    q3 <- (nrow(BNSTDEGs)-BNSTfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    BNSTdf3$OR <- BNSTdf3$pct/BNSTdf3$unifrac
    BNSTdf3$OR.degs <- BNSTdf3$pct.degs/BNSTdf3$unifrac
    BNSTdf4 <- cbind(regname,BNSTdf3,pval)
    BNST.master=rbind(BNST.master, BNSTdf4)

MeA <- data.frame(regulons[i])
    MeA <- as.matrix(MeA)
    regname <- i
    MeADEGsReg <- intersect(MeA, MeADEGs)
    MeADEGsReg <- data.frame(MeADEGsReg)
    Total.MeA = rbind(Total.MeA, MeADEGsReg)
    MeAdect <- intersect(MeA,MeABackground)
    MeAdect <- data.frame(MeAdect)
    MeAfrac <- nrow(MeADEGsReg)
    MeArep <- nrow(MeADEGsReg)/nrow(MeA)
    total <- nrow(MeAdect)
    unifrac <- total/univm
    MeAdf2 <- data.frame(MeAfrac, total)
    MeAdf2$pct <- MeAdf2$MeAfrac/MeAdf2$total
    MeAdf2$pct.degs <- MeAfrac/nrow(MeADEGs)
    MeAdf3 <- cbind(MeAdf2, unifrac)
    q1 <- MeAfrac
    q2 <- (total - MeAfrac)
    q3 <- (nrow(MeADEGs)-MeAfrac)
    q4 <- (univm-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    MeAdf3$OR <- MeAdf3$pct/MeAdf3$unifrac
    MeAdf3$OR.degs <- MeAdf3$pct.degs/MeAdf3$unifrac
    MeAdf4 <- cbind(regname,MeAdf3,pval)
    MeA.master=rbind(MeA.master,MeAdf4)

POA <- data.frame(regulons[i])
    POA <- as.matrix(POA)
    regname <- i
    POADEGsReg <- intersect(POA, POADEGs)
    POADEGsReg <- data.frame(POADEGsReg)
    Total.POA = rbind(Total.POA, POADEGsReg)
    POAdect <- intersect(POA,POABackground)
    POAdect <- data.frame(POAdect)
    POAfrac <- nrow(POADEGsReg)
    POArep <- nrow(POADEGsReg)/nrow(POA)
    total <- nrow(MeAdect)
    unifrac <- total/univp
    POAdf2 <- data.frame(POAfrac, total)
    POAdf2$pct <- POAdf2$POAfrac/POAdf2$total
    POAdf2$pct.degs <- POAfrac/nrow(POADEGs)
    POAdf3 <- cbind(POAdf2, unifrac)
    q1 <- POAfrac
    q2 <- (total - POAfrac)
    q3 <- (nrow(POADEGs)-POAfrac)
    q4 <- (univp-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    POAdf3$OR <- POAdf3$pct/POAdf3$unifrac
    POAdf3$OR.degs <- POAdf3$pct.degs/POAdf3$unifrac
    POAdf4 <- cbind(regname,POAdf3,pval)
    POA.master=rbind(POA.master,POAdf4)

VMH <- data.frame(regulons[i])
    VMH <- as.matrix(VMH)
    regname <- i
    VMHDEGsReg <- intersect(VMH, VMHDEGs)
    VMHDEGsReg <- data.frame(VMHDEGsReg)
    Total.VMH = rbind(Total.VMH, VMHDEGsReg)
    VMHdect <- intersect(VMH,VMHBackground)
    VMHdect <- data.frame(VMHdect)
    VMHfrac <- nrow(VMHDEGsReg)
    VMHrep <- nrow(VMHDEGsReg)/nrow(VMH)
    total <- nrow(VMHdect)
    unifrac <- total/univv
    VMHdf2 <- data.frame(VMHfrac, total)
    VMHdf2$pct <- VMHdf2$VMHfrac/VMHdf2$total
    VMHdf2$pct.degs <- VMHfrac/nrow(VMHDEGs)
    VMHdf3 <- cbind(VMHdf2, unifrac)
    q1 <- VMHfrac
    q2 <- (total - VMHfrac)
    q3 <- (nrow(VMHDEGs)-VMHfrac)
    q4 <- (univv-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    VMHdf3$OR <- VMHdf3$pct/VMHdf3$unifrac
    VMHdf3$OR.degs <- VMHdf3$pct.degs/VMHdf3$unifrac
    VMHdf4 <- cbind(regname,VMHdf3,pval)
    VMH.master=rbind(VMH.master,VMHdf4)

BNSTe <- data.frame(regulons[i])
    BNSTe <- as.matrix(BNSTe)
    regname <- i
    BNSTeDEGsReg <- intersect(BNSTe, BNSTeDEGs)
    BNSTeDEGsReg <- data.frame(BNSTeDEGsReg)
    Total.BNSTe = rbind(Total.BNSTe, BNSTeDEGsReg)
    BNSTedect <- intersect(BNSTe,BNSTBackground)
    BNSTedect <- data.frame(BNSTedect)
    BNSTefrac <- nrow(BNSTeDEGsReg)
    BNSTerep <- nrow(BNSTeDEGsReg)/nrow(BNSTe)
    total <- nrow(BNSTedect)
    unifrac <- total/univb
    BNSTedf2 <- data.frame(BNSTefrac, total)
    BNSTedf2$pct <- BNSTedf2$BNSTefrac/BNSTedf2$total
    BNSTedf2$pct.degs <- BNSTefrac/nrow(BNSTeDEGs)
    BNSTedf3 <- cbind(BNSTedf2, unifrac)
    q1 <- BNSTefrac
    q2 <- (total - BNSTefrac)
    q3 <- (nrow(BNSTeDEGs)-BNSTefrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    BNSTedf3$OR <- BNSTedf3$pct/BNSTedf3$unifrac
    BNSTedf3$OR.degs <- BNSTedf3$pct.degs/BNSTedf3$unifrac
    BNSTedf4 <- cbind(regname,BNSTedf3,pval)
    BNSTe.master=rbind(BNSTe.master,BNSTedf4)
    
MeAe <- data.frame(regulons[i])
    MeAe <- as.matrix(MeAe)
    regname <- i
    MeAeDEGsReg <- intersect(MeAe, MeAeDEGs)
    MeAeDEGsReg <- data.frame(MeAeDEGsReg)
    Total.MeAe = rbind(Total.MeAe, MeAeDEGsReg)
    MeAedect <- intersect(MeAe,MeABackground)
    MeAedect <- data.frame(MeAedect)
    MeAefrac <- nrow(MeAeDEGsReg)
    MeAerep <- nrow(MeAeDEGsReg)/nrow(MeAe)
    total <- nrow(MeAedect)
    unifrac <- total/univm
    MeAedf2 <- data.frame(MeAefrac, total)
    MeAedf2$pct <- MeAedf2$MeAefrac/MeAedf2$total
    MeAedf2$pct.degs <- MeAefrac/nrow(MeAeDEGs)
    MeAedf3 <- cbind(MeAedf2, unifrac)
    q1 <- MeAefrac
    q2 <- (total - MeAefrac)
    q3 <- (nrow(MeAeDEGs)-MeAefrac)
    q4 <- (univm-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    MeAedf3$OR <- MeAedf3$pct/MeAedf3$unifrac
    MeAedf3$OR.degs <- MeAedf3$pct.degs/MeAedf3$unifrac
    MeAedf4 <- cbind(regname,MeAedf3,pval)
    MeAe.master=rbind(MeAe.master,MeAedf4)
    
POAe <- data.frame(regulons[i])
    POAe <- as.matrix(POAe)
    regname <- i
    POAeDEGsReg <- intersect(POAe, POAeDEGs)
    POAeDEGsReg <- data.frame(POAeDEGsReg)
    Total.POAe = rbind(Total.POAe, POAeDEGsReg)
    POAedect <- intersect(POAe,POABackground)
    POAedect <- data.frame(POAedect)
    POAefrac <- nrow(POAeDEGsReg)
    POAerep <- nrow(POAeDEGsReg)/nrow(POAe)
    total <- nrow(POAedect)
    unifrac <- total/univp
    POAedf2 <- data.frame(POAefrac, total)
    POAedf2$pct <- POAedf2$POAefrac/POAedf2$total
    POAedf2$pct.degs <- POAefrac/nrow(POAeDEGs)
    POAedf3 <- cbind(POAedf2, unifrac)
    q1 <- POAefrac
    q2 <- (total - POAefrac)
    q3 <- (nrow(POAeDEGs)-POAefrac)
    q4 <- (univp-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    POAedf3$OR <- POAedf3$pct/POAedf3$unifrac
    POAedf3$OR.degs <- POAedf3$pct.degs/POAedf3$unifrac
    POAedf4 <- cbind(regname,POAedf3,pval)
    POAe.master=rbind(POAe.master,POAedf4)
    
VMHe <- data.frame(regulons[i])
    VMHe <- as.matrix(VMHe)
    regname <- i
    VMHeDEGsReg <- intersect(VMHe, VMHeDEGs)
    VMHeDEGsReg <- data.frame(VMHeDEGsReg)
    Total.VMHe = rbind(Total.VMHe, VMHeDEGsReg)
    VMHedect <- intersect(VMHe,VMHBackground)
    VMHedect <- data.frame(VMHedect)
    VMHefrac <- nrow(VMHeDEGsReg)
    VMHerep <- nrow(VMHeDEGsReg)/nrow(VMHe)
    total <- nrow(VMHedect)
    unifrac <- total/univv
    VMHedf2 <- data.frame(VMHefrac, total)
    VMHedf2$pct <- VMHedf2$VMHefrac/VMHedf2$total
    VMHedf2$pct.degs <- VMHefrac/nrow(VMHeDEGs)
    VMHedf3 <- cbind(VMHedf2, unifrac)
    q1 <- VMHefrac
    q2 <- (total - VMHefrac)
    q3 <- (nrow(VMHeDEGs)-VMHefrac)
    q4 <- (univv-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    VMHedf3$OR <- VMHedf3$pct/VMHedf3$unifrac
    VMHedf3$OR.degs <- VMHedf3$pct.degs/VMHedf3$unifrac
    VMHedf4 <- cbind(regname,VMHedf3,pval)
    VMHe.master=rbind(VMHe.master,VMHedf4)
    
    
Total <- data.frame(regulons[i])
    Total <- as.matrix(Total)
    regname <- i
    TotalDEGsReg <- intersect(Total, DEGs.master)
    TotalDEGsReg <- data.frame(TotalDEGsReg)
    Total.DEGs = rbind(Total.DEGs, TotalDEGsReg)
    Totaldect <- intersect(Total,Background.master)
    Totaldect <- data.frame(Totaldect)
    Totalfrac <- nrow(TotalDEGsReg)
    Toralrep <- nrow(TotalDEGsReg)/nrow(Total)
    total <- nrow(Totaldect)
    unifrac <- total/univt
    Totaldf2 <- data.frame(Totalfrac, total)
    Totaldf2$pct <- Totaldf2$Totalfrac/Totaldf2$total
    Totaldf2$pct.degs <- Totalfrac/nrow(DEGs.master)
    Totaldf3 <- cbind(Totaldf2, unifrac)
    q1 <- Totalfrac
    q2 <- (total - Totalfrac)
    q3 <- (nrow(DEGs.master)-Totalfrac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    Totaldf3$OR <- Totaldf3$pct/Totaldf3$unifrac
    Totaldf3$OR.degs <- Totaldf3$pct.degs/Totaldf3$unifrac
    Totaldf4 <- cbind(regname,Totaldf3,pval)
    Total.master=rbind(Total.master,Totaldf4)

Totale <- data.frame(regulons[i])
    Totale <- as.matrix(Totale)
    regname <- i
    TotaleDEGsReg <- intersect(Total, eDEGs.master)
    TotaleDEGsReg <- data.frame(TotaleDEGsReg)
    Total.eDEGs = rbind(Total.eDEGs, TotaleDEGsReg)
    Totaledect <- intersect(Totale,Background.master)
    Totaledect <- data.frame(Totaledect)
    Totalefrac <- nrow(TotaleDEGsReg)
    Toralerep <- nrow(TotaleDEGsReg)/nrow(Totale)
    total <- nrow(Totaledect)
    unifrac <- total/univt
    Totaledf2 <- data.frame(Totalefrac, total)
    Totaledf2$pct <- Totaledf2$Totalefrac/Totaledf2$total
    Totaledf2$pct.degs <- Totalefrac/nrow(eDEGs.master)
    Totaledf3 <- cbind(Totaledf2, unifrac)
    q1 <- Totalefrac
    q2 <- (total - Totalefrac)
    q3 <- (nrow(eDEGs.master)-Totalefrac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    Totaledf3$OR <- Totaledf3$pct/Totaledf3$unifrac
    Totaledf3$OR.degs <- Totaledf3$pct.degs/Totaledf3$unifrac
    Totaledf4 <- cbind(regname,Totaledf3,pval)
    Totale.master=rbind(Totale.master,Totaledf4)


}
)}

BNSTsig <- BNST.master$pval
B_padj <- p.adjust(BNSTsig, method="BH")
BNST.master <- cbind(BNST.master, B_padj)
write.csv(BNST.master, file=paste0(output,"regulonoverlap_BNSTpregDEGs.csv"))
BNSTgenes <- Total.BNST[,1]
BNSTgenes <- unique(BNSTgenes)
write.csv(BNSTgenes, file=paste0(output,"all_BNST_Degs_In_Regulons.csv"))

MeAsig <- MeA.master$pval
M_padj <- p.adjust(MeAsig, method="BH")
MeA.master <- cbind(MeA.master, M_padj)
write.csv(MeA.master, file=paste0(output,"regulonoverlap_MeApregDEGs.csv"))
MeAgenes <- Total.MeA[,1]
MeAgenes <- unique(MeAgenes)
write.csv(MeAgenes, file=paste0(output,"all_MeA_Degs_In_Regulons.csv"))

POAsig <- POA.master$pval
P_padj <- p.adjust(POAsig, method="BH")
POA.master <- cbind(POA.master, P_padj)
write.csv(POA.master, file=paste0(output,"regulonoverlap_POApregDEGs.csv"))
POAgenes <- Total.POA[,1]
POAgenes <- unique(POAgenes)
write.csv(POAgenes, file=paste0(output,"all_POA_Degs_In_Regulons.csv"))

VMHsig <- VMH.master$pval
V_padj <- p.adjust(VMHsig, method="BH")
VMH.master <- cbind(VMH.master, V_padj)
write.csv(VMH.master, file=paste0(output,"regulonoverlap_VMHpregDEGs.csv"))
VMHgenes <- Total.VMH[,1]
VMHgenes <- unique(VMHgenes)
write.csv(VMHgenes, file=paste0(output,"all_VMH_Degs_In_Regulons.csv"))

Totalsig <- Total.master$pval
T_padj <- p.adjust(Totalsig, method="BH")
Total.master <- cbind(Total.master, T_padj)
write.csv(Total.master, file=paste0(output,"regulonoverlap_AllpregDEGs.csv"))
Totalgenes <- Total.DEGs[,1]
Totalgenes <- unique(Totalgenes)
write.csv(Totalgenes, file=paste0(output,"all_Degs_In_Regulons.csv"))

BNSTesig <- BNSTe.master$pval
Be_padj <- p.adjust(BNSTesig, method="BH")
BNSTe.master <- cbind(BNSTe.master, Be_padj)
write.csv(BNSTe.master, file=paste0(output,"regulonoverlap_BNSTeDEGs.csv"))
BNSTegenes <- Total.BNSTe[,1]
BNSTegenes <- unique(BNSTegenes)
write.csv(BNSTegenes, file=paste0(output,"all_BNST_eDegs_In_Regulons.csv"))

MeAesig <- MeAe.master$pval
Me_padj <- p.adjust(MeAesig, method="BH")
MeAe.master <- cbind(MeAe.master, Me_padj)
write.csv(MeAe.master, file=paste0(output,"regulonoverlap_MeAeDEGs.csv"))
MeAegenes <- Total.MeAe[,1]
MeAegenes <- unique(MeAegenes)
write.csv(MeAegenes, file=paste0(output,"all_MeA_eDegs_In_Regulons.csv"))

POAesig <- POAe.master$pval
Pe_padj <- p.adjust(POAesig, method="BH")
POAe.master <- cbind(POAe.master, Pe_padj)
write.csv(POAe.master, file=paste0(output,"regulonoverlap_POAeDEGs.csv"))
POAegenes <- Total.POAe[,1]
POAegenes <- unique(POAegenes)
write.csv(POAegenes, file=paste0(output,"all_POA_eDegs_In_Regulons.csv"))

VMHesig <- VMHe.master$pval
Ve_padj <- p.adjust(VMHesig, method="BH")
VMHe.master <- cbind(VMHe.master, Ve_padj)
write.csv(VMHe.master, file=paste0(output,"regulonoverlap_VMHeDEGs.csv"))
VMHegenes <- Total.VMHe[,1]
VMHegenes <- unique(VMHegenes)
write.csv(VMHegenes, file=paste0(output,"all_VMH_eDegs_In_Regulons.csv"))

Totalesig <- Totale.master$pval
Te_padj <- p.adjust(Totalesig, method="BH")
Totale.master <- cbind(Totale.master, Te_padj)
write.csv(Totale.master, file=paste0(output,"regulonoverlap_AlleDEGs.csv"))
Totalegenes <- Total.eDEGs[,1]
Totalegenes <- unique(Totalegenes)
write.csv(Totalegenes, file=paste0(output,"all_eDegs_In_Regulons.csv"))
