#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Viper_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viper)
library(tidyverse) 

regdata <- file("ARACHNE/OE_Esr1vst/expressionmtx.txt")
adjfile <- file("ARACHNE/OE_Esr1vst/network_3col.txt")
regul <- aracne2regulon(adjfile, regdata, format="3col", verbose=TRUE)
regul
#write.csv(regul, file="regulon.csv")
mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered.rds")
mySeurat <- SCTransform(mySeurat, vars.to.regres="orig.ident", return.only.var.genes=FALSE)
Male <- subset(x=mySeurat, subset=Hormone=="Intact")
#Male <- subset(x=Male, idents=c("10","16","20","26","7","11","3","21","14","21","4","19","28"))
Primed <- subset(x=mySeurat, subset=Hormone=="Unprimed")
#Primed <- subset(x=Primed, idents=c("10","16","20","26","7","11","3","21","14","21","4","19","28"))
Male.mat <- Male[["SCT"]]@scale.data
Male.mat <- as.matrix(Male.mat)
Primed.mat <- Primed[["SCT"]]@scale.data
Primed.mat <- as.matrix(Primed.mat)
signatures <- viper::rowTtest(x=Male.mat, y=Primed.mat)
#late_activity_matrix_Primed_vs_male <- viperSignature(eset=Primed.mat, ref=Male.mat, repos=TRUE, per=1000, cores=8)
#late_activity_matrix_scaled_Primed_vs_Male <- viper(late_activity_matrix_Primed_vs_male, regulon=regul, cores=8, verbose=F)
#head(late_activity_matrix_scaled_Primed_vs_Male)
#write.csv(late_activity_matrix_scaled_Primed_vs_Male, file="ViperOut/VMH/Vanillavipetest.csv")
signatures <- (qnorm(signatures$p.value/2, lower.tail=FALSE) * sign(signatures$statistic))[,1]
null_model <- viper::ttestNull(x=Male.mat, y=Primed.mat, per=500, repos=TRUE, verbose=T, cores=1)
final_out <- viper::msviper(signatures, regul, null_model, verbose=T, cores=1)
saveRDS(final_out, file="ViperOut/BNST/AllMalevAllPrimedEsr_Inhouse_lomvardas.rds")
df_obj <- summary(final_out, mrs=100)
df_obj
write.csv(df_obj, file="ViperOut/BNST/AllMalevALlPrimed_InhouseLomvardastop100.csv")
pdf(file="ViperOut/BNST/AllMalevAllPrimedInHouseRegmrsplot1Lomvaras.pdf")
plot(final_out, cex=0.7)
dev.off()
