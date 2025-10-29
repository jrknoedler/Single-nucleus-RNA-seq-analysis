#!/usr/bin/env Rscript

library(SCENIC)
library(loomR)
library(AUCell)
library(SCopeLoomR)
library(Matrix)

#Set output directory
output <- "Regulon_Outgroup/VMH_1stpass_Regulons_DefaultNES_Pos"

#Get relevant gene lists and save as matrices
VMHPregup <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_TotalUpPreg.txt")
VMHPregup <- unlist(VMHPregup)
VMHPregup <- as.matrix(VMHPregup)

#Only pregnancy-upregulated genes

VMHDEGs <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_AllvFemale1.5_2.txt")
VMHDEGs <- unlist(VMHDEGs)

VMHBackground <- read.table("Genelists/PregnancyDEGs/UpdatedPregDEGs/VMH_Background.txt")
VMHBackground <- unlist(VMHBackground)


VMHDEGs <- as.matrix(VMHDEGs)

VMHBackground <- as.matrix(VMHBackground)
univv <- nrow(VMHBackground)

path="/scratch/users/knoedler/pySCENIC_Singularity/VMH_Allen/Tests/"

VMH.master = data.frame()

VMH.up.master = data.frame()

Total.VMH = data.frame()

Total.VMH.up = data.frame()

dirs=list.dirs(path, recursive=FALSE)
dirs
for(i in dirs){

try({
loom <- open_loom(paste0(i,"/auc_mtx_filtered.loom"))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)


###total DEGlists
regs <- names(regulons)
for (i in regs){
    
VMH <- data.frame(regulons[i])
    regfolder <- "Regulon_Outgroup/VMH_Regulons/"
    write.csv(VMH, file=paste0(regfolder,i,".csv"))
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
    unifrac <- total/univv
    VMHupdf2 <- data.frame(VMHupfrac, total)
    VMHupdf2$pct <- VMHupdf2$VMHupfrac/VMHupdf2$total
    VMHupdf2$pct.degs <- VMHupfrac/nrow(VMHPregup)
    VMHupdf3 <- cbind(VMHupdf2, unifrac)
    q1 <- VMHupfrac
    q2 <- (total - VMHupfrac)
    q3 <- (nrow(VMHPregup)-VMHupfrac)
    q4 <- (univv-(q3+q2+q1))
    chisq <- cbind(c(q1,q2),c(q3,4))
    t <- chisq.test(chisq)
    pval <- t$p.value
    VMHupdf3$OR <- VMHupdf3$pct/VMHupdf3$unifrac
    VMHupdf3$OR.degs <- VMHupdf3$pct.degs/VMHupdf3$unifrac
    VMHupdf4 <- cbind(regname,VMHupdf3,pval)
    VMH.up.master = rbind(VMH.up.master, VMHupdf4)

loom <- close_loom()
}
})
}

VMHsig <- VMH.master$pval
V_padj <- p.adjust(VMHsig, method="BH")
VMH.master <- cbind(VMH.master, V_padj)
write.csv(VMH.master, file=paste0(output,"regulonoverlap_VMHpregDEGs.csv"))
VMHgenes <- Total.VMH[,1]
VMHgenes <- unique(VMHgenes)
write.csv(VMHgenes, file=paste0(output,"all_VMH_Degs_In_Regulons.csv"))

VMHupsig <- VMH.up.master$pval
VMH_padj <- p.adjust(VMHupsig, method="BH")
VMH.up.master <- cbind(VMH.up.master, VMH_padj)
write.csv(VMH.up.master, file=paste0(output,"regulonoverlap_VMHupDEGs.csv"))
VMHupgenes <- Total.VMH.up[,1]
VMHupgenes <- unique(VMHupgenes)
write.csv(VMHupgenes, file=paste0(output,"all_VMHup_Degs_In_Regulons.csv"))

