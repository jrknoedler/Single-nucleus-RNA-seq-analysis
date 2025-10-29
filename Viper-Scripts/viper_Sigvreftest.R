#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merge_Fig6storyboard_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viper)
library(tidyverse) 

regdata <- file("ARACHNE/AllBulkvst/Expressionmtx.txt")
adjfile <- file("ARACHNE/AllBulkvst/network_3col.txt")
regul <- aracne2regulon(adjfile, regdata, format="3col", verbose=TRUE)
mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered3_Paperanalysis_30pcs.rds")

#BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
Male <- subset(x=mySeurat, subset=Hormone=="Intact")
#Core <- subset(x=mySeurat, idents=c("0","15","16","3","4","8","5","13","28","24","25","22","27","6"))
Primed <- subset(x=mySeurat, subset=Hormone=="Primed")
#Shell <- subset(x=mySeurat, idents=c("19","11","2","17","7","14","26","19"))
Ref <- subset(x=mySeurat, idents=c("0"), invert=TRUE)
#Male.mat <- Male[["SCT"]]@scale.data
Primed0 <- subset(x=Primed, idents=c("0"))
Male0 <- subset(x=Male, idents=c("0"))
Primed0.mat <- Primed0$SCT@scale.data
Male0.mat <- Male0$SCT@scale.data
Ref.mat <- Ref$SCT@scale.data
MaleSig <- viperSignature(eset=Male0.mat, ref=Ref.mat, per=1000, cores=8)
PrimedSig <- viperSignature(eset=Primed0.mat, ref=Ref.mat, per=1000, cores=8)
Male_act <- viper(MaleSig, regulon=regul, cores=8 ,verbose=T)
Primed_act <- viper(PrimedSig, regulon=regul, cores=8, verbose=T)
Male_act
Primed_act
#Male.mat <- as.matrix(Male.mat)
#Primed.mat <- Primed[["SCT"]]@scale.data
#Primed.mat <- as.matrix(Primed.mat)
#Core.mat <- Core$SCT@scale.data
#Shell.mat <- Shell$SCT@scale.data
signatures <- viper::rowTtest(x=Core.mat, y=Shell.mat)
#late_activity_matrix_Primed_vs_male <- viperSignature(eset=Primed.mat, ref=Male.mat, repos=TRUE, per=1000, cores=8)
#late_activity_matrix_scaled_Primed_vs_Male <- viper(late_activity_matrix_Primed_vs_male, regulon=regul, cores=8, verbose=F)
#head(late_activity_matrix_scaled_Primed_vs_Male)
#write.csv(late_activity_matrix_scaled_Primed_vs_Male, file="ViperOut/VMH/Vanillavipetest.csv")
#signatures <- (qnorm(signatures$p.value/2, lower.tail=FALSE) * sign(signatures$statistic))[,1]
#null_model <- viper::ttestNull(x=Core.mat, y=Shell.mat, per=1000, repos=TRUE, verbose=T, cores=8)
#final_out <- viper::msviper(signatures, regul, null_model, verbose=T, cores=8)
#saveRDS(final_out, file="ViperOut/VMH/Haniscripttest_1000permutecorrectcores.rds")
#df_obj <- summary(final_out, length=(final_out$es$nes)) %>% as.data.frame() %>% arrange(FDR)
#final_out
#df_obj