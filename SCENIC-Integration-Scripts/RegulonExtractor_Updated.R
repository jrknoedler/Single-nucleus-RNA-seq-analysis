#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(AUCell)
library(SCopeLoomR)
library(Matrix)

#Set output directory
output <- "Seurat/Pregnancy_Full_Dataset_Analysis/Regulons_1stpass/500bp_"

#Get relevant gene lists and save as matrices
BNSTDEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_AllvFemale1.5_3.txt")
BNSTDEGs <- unlist(BNSTDEGs)
BNSTDEGs <- as.matrix(BNSTDEGs)
BNSTBackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_Background.txt")
BNSTBackground <- as.matrix(BNSTBackground)
univb <- nrow(BNSTBackground)

MeADEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_AllvFemale1.5_2.txt")
MeADEGs <- unlist(MeADEGs)
MeADEGs <- as.matrix(MeADEGs)
MeABackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_Background.txt")
MeABackground <- as.matrix(MeABackground)
univm <- nrow(MeABackground)

POADEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_AllvFemale1.5_2.txt")
POADEGs <- unlist(POADEGs)
POADEGs <- as.matrix(POADEGs)
POABackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_Background.txt")
POABackground <- as.matrix(POABackground)
univp <- nrow(POABackground)

VMHDEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_AllvFemale1.5_2.txt")
VMHDEGs <- unlist(VMHDEGs)
VMHDEGs <- as.matrix(VMHDEGs)
VMHBackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_Background.txt")
VMHBackground <- as.matrix(VMHBackground)
univv <- nrow(VMHBackground)

path="/scratch/users/knoedler/pySCENIC_Singularity/Subsampled_TFs_OnlyPos_500bp"
BNST.master =data.frame()
MeA.master = data.frame()
POA.master = data.frame()
VMH.master = data.frame()
dirs=list.dirs(path, recursive=FALSE)
dirs
for(i in dirs){

try({
loom <- open_loom(paste0(i,"/auc_mtx_filtered.loom"))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

regs <- names(regulons)
for (i in regs){
    BNST <- data.frame(regulons[i])
    regfolder <- "Seurat/Pregnancy_Full_Dataset_Analysis/Regulons_1stpass/Regulons_500bpOnly/"
    write.csv(BNST, file=paste0(regfolder,i,".csv"))
    BNST <- as.matrix(BNST)
    regname <- i
    BNSTDEGsReg <- intersect(BNST, BNSTDEGs)
    BNSTDEGsReg <- data.frame(BNSTDEGsReg)
    BNSTdect <- intersect(BNST,BNSTBackground)
    BNSTdect <- data.frame(BNSTdect)
    BNSTfrac <- nrow(BNSTDEGsReg)
    total <- nrow(BNSTdect)
    unifrac <- total/univb
    BNSTdf2 <- data.frame(BNSTfrac, total)
    BNSTdf2$pct <- BNSTdf2$BNSTfrac/BNSTdf2$total
    BNSTdf3 <- cbind(BNSTdf2, unifrac)
    q1 <- BNSTfrac
    q2 <- (total - BNSTfrac)
    q3 <- (nrow(BNSTDEGs)-BNSTfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    BNSTdf3$OR <- BNSTdf3$pct/BNSTdf3$unifrac
    BNSTdf4 <- cbind(regname,BNSTdf3,pval)
    BNST.master=rbind(BNST.master, BNSTdf4)

MeA <- data.frame(regulons[i])
    MeA <- as.matrix(MeA)
    regname <- i
    MeADEGsReg <- intersect(MeA, MeADEGs)
    MeADEGsReg <- data.frame(MeADEGsReg)
    MeAdect <- intersect(MeA,MeABackground)
    MeAdect <- data.frame(MeAdect)
    MeAfrac <- nrow(MeADEGsReg)
    total <- nrow(MeAdect)
    unifrac <- total/univm
    MeAdf2 <- data.frame(MeAfrac, total)
    MeAdf2$pct <- MeAdf2$MeAfrac/MeAdf2$total
    MeAdf3 <- cbind(MeAdf2, unifrac)
    q1 <- MeAfrac
    q2 <- (total - MeAfrac)
    q3 <- (nrow(MeADEGs)-MeAfrac)
    q4 <- (univm-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    MeAdf3$OR <- MeAdf3$pct/MeAdf3$unifrac
    MeAdf4 <- cbind(regname,MeAdf3,pval)
    MeA.master=rbind(MeA.master,MeAdf4)

POA <- data.frame(regulons[i])
    POA <- as.matrix(POA)
    regname <- i
    POADEGsReg <- intersect(POA, POADEGs)
    POADEGsReg <- data.frame(POADEGsReg)
    POAdect <- intersect(POA,POABackground)
    POAdect <- data.frame(POAdect)
    POAfrac <- nrow(POADEGsReg)
    total <- nrow(MeAdect)
    unifrac <- total/univp
    POAdf2 <- data.frame(POAfrac, total)
    POAdf2$pct <- POAdf2$POAfrac/POAdf2$total
    POAdf3 <- cbind(POAdf2, unifrac)
    q1 <- POAfrac
    q2 <- (total - POAfrac)
    q3 <- (nrow(POADEGs)-POAfrac)
    q4 <- (univp-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    POAdf3$OR <- POAdf3$pct/POAdf3$unifrac
    POAdf4 <- cbind(regname,POAdf3,pval)
    POA.master=rbind(POA.master,POAdf4)

VMH <- data.frame(regulons[i])
    VMH <- as.matrix(VMH)
    regname <- i
    VMHDEGsReg <- intersect(VMH, VMHDEGs)
    VMHDEGsReg <- data.frame(VMHDEGsReg)
    VMHdect <- intersect(VMH,VMHBackground)
    VMHdect <- data.frame(VMHdect)
    VMHfrac <- nrow(VMHDEGsReg)
    total <- nrow(VMHdect)
    unifrac <- total/univv
    VMHdf2 <- data.frame(VMHfrac, total)
    VMHdf2$pct <- VMHdf2$VMHfrac/VMHdf2$total
    VMHdf3 <- cbind(VMHdf2, unifrac)
    q1 <- VMHfrac
    q2 <- (total - VMHfrac)
    q3 <- (nrow(VMHDEGs)-VMHfrac)
    q4 <- (univv-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    VMHdf3$OR <- VMHdf3$pct/VMHdf3$unifrac
    VMHdf4 <- cbind(regname,VMHdf3,pval)
    VMH.master=rbind(VMH.master,VMHdf4)

loom <- close_loom()
}
})
}

BNSTsig <- BNST.master$pval
B_padj <- p.adjust(BNSTsig, method="BH")
BNST.master <- cbind(BNST.master, B_padj)
write.csv(BNST.master, file=paste0(output,"regulonoverlap_BNSTpregDEGs.csv"))

MeAsig <- MeA.master$pval
M_padj <- p.adjust(MeAsig, method="BH")
MeA.master <- cbind(MeA.master, M_padj)
write.csv(MeA.master, file=paste0(output,"regulonoverlap_MeApregDEGs.csv"))

POAsig <- POA.master$pval
P_padj <- p.adjust(POAsig, method="BH")
POA.master <- cbind(POA.master, P_padj)
write.csv(POA.master, file=paste0(output,"regulonoverlap_POApregDEGs.csv"))

VMHsig <- VMH.master$pval
V_padj <- p.adjust(VMHsig, method="BH")
VMH.master <- cbind(VMH.master, V_padj)
write.csv(VMH.master, file=paste0(output,"regulonoverlap_VMHpregDEGs.csv"))

