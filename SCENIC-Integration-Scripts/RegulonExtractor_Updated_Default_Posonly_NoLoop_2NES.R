#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(AUCell)
library(SCopeLoomR)
library(Matrix)

#Set output directory
output <- "Seurat/Pregnancy_Full_Dataset_Analysis/Regulons_Pseudobulk/2NES_Onlypos/"

#Get relevant gene lists and save as matrices
#Raw DEG lists
BNSTPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_TotalUpPreg.txt")
BNSTPregup <- unlist(BNSTPregup)
BNSTPregup <- as.matrix(BNSTPregup)
MeAPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_TotalUpPreg.txt")
MeAPregup <- unlist(MeAPregup)
MeAPregup <- as.matrix(MeAPregup)
POAPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_TotalUpPreg.txt")
POAPregup <- unlist(POAPregup)
POAPregup <- as.matrix(POAPregup)
VMHPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_TotalUpPreg.txt")
VMHPregup <- unlist(VMHPregup)
VMHPregup <- as.matrix(VMHPregup)

#Only pregnancy-upregulated genes
BNSTDEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_AllvFemale1.5_3.txt")
BNSTDEGs <- unlist(BNSTDEGs)
MeADEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_AllvFemale1.5_2.txt")
MeADEGs <- unlist(MeADEGs)
POADEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_AllvFemale1.5_2.txt")
POADEGs <- unlist(POADEGs)
VMHDEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_AllvFemale1.5_2.txt")
VMHDEGs <- unlist(VMHDEGs)
DEGs.master1 <- union(BNSTDEGs,MeADEGs)
DEGs.master2 <- union(DEGs.master1, POADEGs)
DEGs.master3 <- union(DEGs.master2, VMHDEGs)
DEGs.master <- as.matrix(DEGs.master3)

#Unified K clusters
K1.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K1_all.txt")
K1.DEGs <- unlist(K1.DEGs)
K1.DEGs <- as.matrix(K1.DEGs)
K2.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K2_all.txt")
K2.DEGs <- unlist(K2.DEGs)
K2.DEGs <- as.matrix(K2.DEGs)
K3.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K3_all.txt")
K3.DEGs <- unlist(K3.DEGs)
K3.DEGs <- as.matrix(K3.DEGs)
head(K3.DEGs)
dim(K3.DEGs)
K4.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K4_all.txt")
K4.DEGs <- unlist(K4.DEGs)
K4.DEGs <- as.matrix(K4.DEGs)
K5.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K5_all.txt")
K5.DEGs <- unlist(K5.DEGs)
K5.DEGs <- as.matrix(K5.DEGs)
K6.DEGs <- read.table("Kmeans_Reassigned/Collapsed_Kclusters/Female_K6_all.txt")
K6.DEGs <- unlist(K6.DEGs)
K6.DEGs <- as.matrix(K6.DEGs)


BNSTBackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/BNST_Background.txt")
BNSTBackground <- unlist(BNSTBackground)
MeABackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/MeA_Background.txt")
MeABackground <- unlist(MeABackground)
POABackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/POA_Background.txt")
POABackground <- unlist(POABackground)
VMHBackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_Background.txt")
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

path="/scratch/users/knoedler/pySCENIC_Singularity/Subsampled_TFs_Fullmatrixtest"
BNST.master =data.frame()
MeA.master = data.frame()
POA.master = data.frame()
VMH.master = data.frame()
Total.master = data.frame()
K1.master = data.frame()
K2.master = data.frame()
K3.master = data.frame()
K4.master = data.frame()
K5.master = data.frame()
K6.master = data.frame()
BNST.up.master = data.frame()
MeA.up.master = data.frame()
POA.up.master = data.frame()
VMH.up.master = data.frame()
Total.BNST = data.frame()
Total.MeA = data.frame()
Total.POA = data.frame()
Total.VMH = data.frame()
Total.DEGs = data.frame()
Total.BNST.up = data.frame()
Total.MeA.up = data.frame()
Total.POA.up = data.frame()
Total.VMH.up = data.frame()
Total.K1 = data.frame()
Total.K2 = data.frame()
Total.K3 = data.frame()
Total.K4 = data.frame()
Total.K5 = data.frame()
Total.K6 = data.frame()




loom <- open_loom("pySCENIC_Singularity/PseudobulkGRN/GRN_Output/auc_mtx_filtered_2NES_onlyPos_Pseudobulkregulons.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)


###total DEGlists
regs <- names(regulons)
for (i in regs){
try({
    BNST <- data.frame(regulons[i])
    regfolder <- "Seurat/Pregnancy_Full_Dataset_Analysis/Regulons_Pseudobulk/2NES_Onlypos/Regulons/"
    write.csv(BNST, file=paste0(regfolder,i,".csv"))
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

###Pregnancy-upregulated genes
BNSTup <- data.frame(regulons[i])
    BNSTup <- as.matrix(BNSTup)
    regname <- i
    BNSTupDEGsReg <- intersect(BNSTup, BNSTPregup)
    BNSTupDEGsReg <- data.frame(BNSTupDEGsReg)
    Total.BNST.up = rbind(Total.BNST.up, BNSTupDEGsReg)
    BNSTupdect <- intersect(BNSTup,BNSTBackground)
    BNSTupdect <- data.frame(BNSTupdect)
    BNSTupfrac <- nrow(BNSTupDEGsReg)
    BNSTuprep <- nrow(BNSTupDEGsReg)/nrow(BNSTup)
    total <- nrow(BNSTupdect)
    unifrac <- total/univb
    BNSTupdf2 <- data.frame(BNSTupfrac, total)
    BNSTupdf2$pct <- BNSTupdf2$BNSTupfrac/BNSTupdf2$total
    BNSTupdf2$pct.degs <- BNSTupfrac/nrow(BNSTPregup)
    BNSTupdf3 <- cbind(BNSTupdf2, unifrac)
    q1 <- BNSTupfrac
    q2 <- (total - BNSTupfrac)
    q3 <- (nrow(BNSTPregup)-BNSTupfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    BNSTupdf3$OR <- BNSTupdf3$pct/BNSTupdf3$unifrac
    BNSTupdf3$OR.degs <- BNSTupdf3$pct.degs/BNSTupdf3$unifrac
    BNSTupdf4 <- cbind(regname,BNSTupdf3,pval)
    BNST.up.master = rbind(BNST.up.master, BNSTupdf4)

MeAup <- data.frame(regulons[i])
    MeAup <- as.matrix(MeAup)
    regname <- i
    MeAupDEGsReg <- intersect(MeAup, MeAPregup)
    MeAupDEGsReg <- data.frame(MeAupDEGsReg)
    Total.MeA.up = rbind(Total.MeA.up, MeAupDEGsReg)
    MeAupdect <- intersect(MeAup,MeABackground)
    MeAupdect <- data.frame(MeAupdect)
    MeAupfrac <- nrow(MeAupDEGsReg)
    MeAuprep <- nrow(MeAupDEGsReg)/nrow(MeAup)
    total <- nrow(MeAupdect)
    unifrac <- total/univb
    MeAupdf2 <- data.frame(MeAupfrac, total)
    MeAupdf2$pct <- MeAupdf2$MeAupfrac/MeAupdf2$total
    MeAupdf2$pct.degs <- MeAupfrac/nrow(MeAPregup)
    MeAupdf3 <- cbind(MeAupdf2, unifrac)
    q1 <- MeAupfrac
    q2 <- (total - MeAupfrac)
    q3 <- (nrow(MeAPregup)-MeAupfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    MeAupdf3$OR <- MeAupdf3$pct/MeAupdf3$unifrac
    MeAupdf3$OR.degs <- MeAupdf3$pct.degs/MeAupdf3$unifrac
    MeAupdf4 <- cbind(regname,MeAupdf3,pval)
    MeA.up.master = rbind(MeA.up.master, MeAupdf4)

POAup <- data.frame(regulons[i])
    POAup <- as.matrix(POAup)
    regname <- i
    POAupDEGsReg <- intersect(POAup, POAPregup)
    POAupDEGsReg <- data.frame(POAupDEGsReg)
    Total.POA.up = rbind(Total.POA.up, POAupDEGsReg)
    POAupdect <- intersect(POAup,POABackground)
    POAupdect <- data.frame(POAupdect)
    POAupfrac <- nrow(POAupDEGsReg)
    POAuprep <- nrow(POAupDEGsReg)/nrow(POAup)
    total <- nrow(POAupdect)
    unifrac <- total/univb
    POAupdf2 <- data.frame(POAupfrac, total)
    POAupdf2$pct <- POAupdf2$POAupfrac/POAupdf2$total
    POAupdf2$pct.degs <- POAupfrac/nrow(POAPregup)
    POAupdf3 <- cbind(POAupdf2, unifrac)
    q1 <- POAupfrac
    q2 <- (total - POAupfrac)
    q3 <- (nrow(POAPregup)-POAupfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    POAupdf3$OR <- POAupdf3$pct/POAupdf3$unifrac
    POAupdf3$OR.degs <- POAupdf3$pct.degs/POAupdf3$unifrac
    POAupdf4 <- cbind(regname,POAupdf3,pval)
    POA.up.master = rbind(POA.up.master, POAupdf4)

VMHup <- data.frame(regulons[i])
    VMHup <- as.matrix(VMHup)
    regname <- i
    VMHupDEGsReg <- intersect(VMHup, VMHPregup)
    VMHupDEGsReg <- data.frame(VMHupDEGsReg)
    Total.VMH.up = rbind(Total.VMH.up, VMHupDEGsReg)
    VMHupdect <- intersect(VMHup,VMHBackground)
    VMHupdect <- data.frame(VMHupdect)
    VMHupfrac <- nrow(VMHupDEGsReg)
    VMHuprep <- nrow(VMHupDEGsReg)/nrow(VMHup)
    total <- nrow(VMHupdect)
    unifrac <- total/univb
    VMHupdf2 <- data.frame(VMHupfrac, total)
    VMHupdf2$pct <- VMHupdf2$VMHupfrac/VMHupdf2$total
    VMHupdf2$pct.degs <- VMHupfrac/nrow(VMHPregup)
    VMHupdf3 <- cbind(VMHupdf2, unifrac)
    q1 <- VMHupfrac
    q2 <- (total - VMHupfrac)
    q3 <- (nrow(VMHPregup)-VMHupfrac)
    q4 <- (univb-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    VMHupdf3$OR <- VMHupdf3$pct/VMHupdf3$unifrac
    VMHupdf3$OR.degs <- VMHupdf3$pct.degs/VMHupdf3$unifrac
    VMHupdf4 <- cbind(regname,VMHupdf3,pval)
    VMH.up.master = rbind(VMH.up.master, VMHupdf4)

###Kmeans lists
K1 <- data.frame(regulons[i])
    K1 <- as.matrix(K1)
    regname <- i
    K1DEGsReg <- intersect(K1, K1.DEGs)
    K1DEGsReg <- data.frame(K1DEGsReg)
    Total.K1 = rbind(Total.K1, K1DEGsReg)
    K1dect <- intersect(K1,Background.master)
    K1dect <- data.frame(K1dect)
    K1frac <- nrow(K1DEGsReg)
    K1rep <- nrow(K1DEGsReg)/nrow(K1)
    total <- nrow(K1dect)
    unifrac <- total/univt
    K1df2 <- data.frame(K1frac, total)
    K1df2$pct <- K1df2$K1frac/K1df2$total
    K1df2$pct.degs <- K1frac/nrow(K1.DEGs)
    K1df3 <- cbind(K1df2, unifrac)
    q1 <- K1frac
    q2 <- (total - K1frac)
    q3 <- (nrow(K1.DEGs)-K1frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K1df3$OR <- K1df3$pct/K1df3$unifrac
    K1df3$OR.degs <- K1df3$pct.degs/K1df3$unifrac
    K1df4 <- cbind(regname,K1df3,pval)
    K1.master=rbind(K1.master,K1df4)

K2 <- data.frame(regulons[i])
    K2 <- as.matrix(K2)
    regname <- i
    K2DEGsReg <- intersect(K2,K2.DEGs)
    K2DEGsReg <- data.frame(K2DEGsReg)
    Total.K2 = rbind(Total.K2, K2DEGsReg)
    K2dect <- intersect(K2,Background.master)
    K2dect <- data.frame(K2dect)
    K2frac <- nrow(K2DEGsReg)
    K2rep <- nrow(K1DEGsReg)/nrow(K2)
    total <- nrow(K2dect)
    unifrac <- total/univt
    K2df2 <- data.frame(K2frac, total)
    K2df2$pct <- K2df2$K2frac/K2df2$total
    K2df2$pct.degs <- K2frac/nrow(K2.DEGs)
    K2df3 <- cbind(K2df2, unifrac)
    q1 <- K2frac
    q2 <- (total - K2frac)
    q3 <- (nrow(K2.DEGs)-K2frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K2df3$OR <- K2df3$pct/K2df3$unifrac
    K2df3$OR.degs <- K2df3$pct.degs/K2df3$unifrac
    K2df4 <- cbind(regname,K2df3,pval)
    K2.master=rbind(K2.master,K2df4)

K3 <- data.frame(regulons[i])
    K3 <- as.matrix(K3)
    regname <- i
    K3DEGsReg <- intersect(K3, K3.DEGs)
    K3DEGsReg <- data.frame(K3DEGsReg)
    Total.K3 = rbind(Total.K3, K3DEGsReg)
    K3dect <- intersect(K3,Background.master)
    K3dect <- data.frame(K3dect)
    K3frac <- nrow(K3DEGsReg)
    K3rep <- nrow(K3DEGsReg)/nrow(K3)
    total <- nrow(K3dect)
    unifrac <- total/univt
    K3df2 <- data.frame(K3frac, total)
    K3df2$pct <- K3df2$K3frac/K3df2$total
    K3df2$pct.degs <- K3frac/nrow(K3.DEGs)
    K3df3 <- cbind(K3df2, unifrac)
    q1 <- K3frac
    q2 <- (total - K3frac)
    q3 <- (nrow(K3.DEGs)-K3frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K3df3$OR <- K3df3$pct/K3df3$unifrac
    K3df3$OR.degs <- K3df3$pct.degs/K3df3$unifrac
    K3df4 <- cbind(regname,K3df3,pval)
    K3.master=rbind(K3.master,K3df4)

K4 <- data.frame(regulons[i])
    K4 <- as.matrix(K4)
    regname <- i
    K4DEGsReg <- intersect(K4,K4.DEGs)
    K4DEGsReg <- data.frame(K4DEGsReg)
    Total.K4 = rbind(Total.K4, K4DEGsReg)
    K4dect <- intersect(K4,Background.master)
    K4dect <- data.frame(K4dect)
    K4frac <- nrow(K4DEGsReg)
    K4rep <- nrow(K4DEGsReg)/nrow(K4)
    total <- nrow(K4dect)
    unifrac <- total/univt
    K4df2 <- data.frame(K4frac, total)
    K4df2$pct <- K4df2$K4frac/K4df2$total
    K4df2$pct.degs <- K4frac/nrow(K4.DEGs)
    K4df3 <- cbind(K4df2, unifrac)
    q1 <- K4frac
    q2 <- (total - K4frac)
    q3 <- (nrow(K4.DEGs)-K4frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K4df3$OR <- K4df3$pct/K4df3$unifrac
    K4df3$OR.degs <- K4df3$pct.degs/K4df3$unifrac
    K4df4 <- cbind(regname,K4df3,pval)
    K4.master=rbind(K4.master,K4df4)

K5 <- data.frame(regulons[i])
    K5 <- as.matrix(K5)
    regname <- i
    K5DEGsReg <- intersect(K5,K5.DEGs)
    K5DEGsReg <- data.frame(K5DEGsReg)
    Total.K5 = rbind(Total.K5, K5DEGsReg)
    K5dect <- intersect(K5,Background.master)
    K5dect <- data.frame(K5dect)
    K5frac <- nrow(K5DEGsReg)
    K5rep <- nrow(K5DEGsReg)/nrow(K5)
    total <- nrow(K5dect)
    unifrac <- total/univt
    K5df2 <- data.frame(K5frac, total)
    K5df2$pct <- K5df2$K5frac/K5df2$total
    K5df2$pct.degs <- K5frac/nrow(K5.DEGs)
    K5df3 <- cbind(K5df2, unifrac)
    q1 <- K5frac
    q2 <- (total - K5frac)
    q3 <- (nrow(K5.DEGs)-K5frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K5df3$OR <- K5df3$pct/K5df3$unifrac
    K5df3$OR.degs <- K5df3$pct.degs/K5df3$unifrac
    K5df4 <- cbind(regname,K5df3,pval)
    K5.master=rbind(K5.master,K5df4)

K6 <- data.frame(regulons[i])
    K6 <- as.matrix(K6)
    regname <- i
    K6DEGsReg <- intersect(K6, K6.DEGs)
    K6DEGsReg <- data.frame(K6DEGsReg)
    Total.K6 = rbind(Total.K6, K6DEGsReg)
    K6dect <- intersect(K6,Background.master)
    K6dect <- data.frame(K6dect)
    K6frac <- nrow(K6DEGsReg)
    K6rep <- nrow(K6DEGsReg)/nrow(K6)
    total <- nrow(K6dect)
    unifrac <- total/univt
    K6df2 <- data.frame(K6frac, total)
    K6df2$pct <- K6df2$K6frac/K6df2$total
    K6df2$pct.degs <- K6frac/nrow(K6.DEGs)
    K6df3 <- cbind(K6df2, unifrac)
    q1 <- K6frac
    q2 <- (total - K6frac)
    q3 <- (nrow(K6.DEGs)-K6frac)
    q4 <- (univt-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    K6df3$OR <- K6df3$pct/K6df3$unifrac
    K6df3$OR.degs <- K6df3$pct.degs/K6df3$unifrac
    K6df4 <- cbind(regname,K6df3,pval)
    K6.master=rbind(K6.master,K6df4)


loom <- close_loom()
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

BNSTupsig <- BNST.up.master$pval
BP_padj <- p.adjust(BNSTupsig, method="BH")
BNST.up.master <- cbind(BNST.up.master, BP_padj)
write.csv(BNST.up.master, file=paste0(output,"regulonoverlap_BNSTupDEGs.csv"))
BNSTupgenes <- Total.BNST.up[,1]
BNSTupgenes <- unique(BNSTupgenes)
write.csv(BNSTupgenes, file=paste0(output,"all_BNSTup_Degs_In_Regulons.csv"))

MeAupsig <- MeA.up.master$pval
MP_padj <- p.adjust(MeAupsig, method="BH")
MeA.up.master <- cbind(MeA.up.master, MP_padj)
write.csv(MeA.up.master, file=paste0(output,"regulonoverlap_MeAupDEGs.csv"))
MeAupgenes <- Total.MeA.up[,1]
MeAupgenes <- unique(MeAupgenes)
write.csv(MeAupgenes, file=paste0(output,"all_MeAup_Degs_In_Regulons.csv"))

POAupsig <- POA.up.master$pval
POA_padj <- p.adjust(POAupsig, method="BH")
POA.up.master <- cbind(POA.up.master, POA_padj)
write.csv(POA.up.master, file=paste0(output,"regulonoverlap_POAupDEGs.csv"))
POAupgenes <- Total.POA.up[,1]
POAupgenes <- unique(POAupgenes)
write.csv(POAupgenes, file=paste0(output,"all_POAup_Degs_In_Regulons.csv"))

VMHupsig <- VMH.up.master$pval
VMH_padj <- p.adjust(VMHupsig, method="BH")
VMH.up.master <- cbind(VMH.up.master, BP_padj)
write.csv(VMH.up.master, file=paste0(output,"regulonoverlap_VMHupDEGs.csv"))
VMHupgenes <- Total.VMH.up[,1]
VMHupgenes <- unique(VMHupgenes)
write.csv(VMHupgenes, file=paste0(output,"all_VMHup_Degs_In_Regulons.csv"))

Totalsig <- Total.master$pval
T_padj <- p.adjust(Totalsig, method="BH")
Total.master <- cbind(Total.master, T_padj)
write.csv(Total.master, file=paste0(output,"regulonoverlap_AllpregDEGs.csv"))
Totalgenes <- Total.DEGs[,1]
Totalgenes <- unique(Totalgenes)
write.csv(Totalgenes, file=paste0(output,"all_Degs_In_Regulons.csv"))

K1sig <- K1.master$pval
K1_padj <- p.adjust(K1sig, method="BH")
K1.master <- cbind(K1.master, K1_padj)
write.csv(K1.master, file=paste0(output,"regulonoverlap_AllK1.csv"))
K1genes <- Total.K1[,1]
K1genes <- unique(K1genes)
write.csv(K1genes, file=paste0(output,"all_K1_Degs_In_Regulons.csv"))

K2sig <- K2.master$pval
K2_padj <- p.adjust(K2sig, method="BH")
K2.master <- cbind(K2.master, K2_padj)
write.csv(K2.master, file=paste0(output,"regulonoverlap_AllK2.csv"))
K2genes <- Total.K2[,1]
K2genes <- unique(K2genes)
write.csv(K2genes, file=paste0(output,"all_K2_Degs_In_Regulons.csv"))

K3sig <- K3.master$pval
K3_padj <- p.adjust(K3sig, method="BH")
K3.master <- cbind(K3.master, K3_padj)
write.csv(K3.master, file=paste0(output,"regulonoverlap_AllK3.csv"))
K3genes <- Total.K3[,1]
K3genes <- unique(K3genes)
write.csv(K3genes, file=paste0(output,"all_K3_Degs_In_Regulons.csv"))

K4sig <- K4.master$pval
K4_padj <- p.adjust(K4sig, method="BH")
K4.master <- cbind(K4.master, K4_padj)
write.csv(K4.master, file=paste0(output,"regulonoverlap_AllK4.csv"))
head(Total.K4)
K4genes <- Total.K4[,1]
K4genes <- unique(K4genes)
write.csv(K4genes, file=paste0(output,"all_K4_Degs_In_Regulons.csv"))

K5sig <- K5.master$pval
K5_padj <- p.adjust(K5sig, method="BH")
K5.master <- cbind(K5.master, K5_padj)
write.csv(K5.master, file=paste0(output,"regulonoverlap_AllK5.csv"))
K5genes <- Total.K5[,1]
K5genes <- unique(K5genes)
write.csv(K5genes, file=paste0(output,"all_K5_Degs_In_Regulons.csv"))

K6sig <- K6.master$pval
K6_padj <- p.adjust(K6sig, method="BH")
K6.master <- cbind(K6.master, K6_padj)
write.csv(K6.master, file=paste0(output,"regulonoverlap_AllK6.csv"))
K6genes <- Total.K6[,1]
K6genes <- unique(K6genes)
write.csv(K6genes, file=paste0(output,"all_K6_Degs_In_Regulons.csv"))