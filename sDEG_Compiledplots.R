#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
data <- read.table("AdultVMH_MvP_Tpm_Log10pluseone.txt", header=TRUE)

VMHMvP <- ggplot(data, aes(x=Male, y=Primed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","deeppink1,"Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "black", color = "black", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="black"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="white"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="white", size=20)) + theme(axis.title.y = element_text(color="white", size=20))+theme(axis.text.x=element_text(size=16, color="white"), axis.text.y = element_text(color="white", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)+ geom_text_repel(data=subset(data, Gene=="Cckar"), label="Cckar", color="white", nudge_x=-0.5, nudge_y=0.5, size=16)

data <- read.table("BNST_MvU_1.5tpmDraft.txt", header=TRUE)


BNSTMvU <- ggplot(data, aes(x=Male, y=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5)

data <- read.table("MeA_MvU_draft1.5.txt", header=TRUE)

MeAMvU <- ggplot(data, aes(x=Male, y=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) +xlim(0,5)+ylim(0,5)

data <- read.table("POA_MvU_draft1.5.txt", header=TRUE)

POAMvU <- ggplot(data, aes(x=Male, y=Unprimed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)

data <- read.table("VMH_MvU_1.5draft.txt", header=TRUE)

VMHMvU <- ggplot(data, aes(x=Male, y=Diestrus)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("Blue","Gray","Green"))+scale_alpha_manual(values=c(1,0.1,1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Log2 Counts (Male)", y="Log2 Counts (Unprimed Female)") + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)

data <- read.table("BNST_MvP_TPM_Log10pluseone.txt", header=TRUE)

BNSTMvP <- ggplot(data, aes(x=Male, y=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) 

data <- read.table("MeA_MvP_again.txt", header=TRUE)

MeAMvP <- ggplot(data, aes(x=Male, y=Primed)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1) + ylim(min=0,5) 

data <- read.table("AdultPOA_MvP_Tpm.csv", header=TRUE)

POAMvP <- ggplot(data, aes(x=Male, y=Female)) + geom_point(aes(color=Category,alpha=Category))+scale_color_manual(values=c("deeppink1","Blue","Gray"))+scale_alpha_manual(values=c(1,1,0.1))+theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+theme(plot.background= element_rect(fill="white"))+theme(axis.line=element_line(size=0.5, linetype="solid", color="black"))+labs(x="Male TPM (Log10+1)", y="Female TPM (Log10+1)") + xlim(0,5)+ylim(0,5) + theme(axis.title.x = element_text(color="black", size=20)) + theme(axis.title.y = element_text(color="black", size=20))+theme(axis.text.x=element_text(size=16, color="black"), axis.text.y = element_text(color="black", size=16))+theme(legend.position="none")+ coord_equal(ratio=1)

pdf(file=paste0("sDEG_Fig1compiled.pdf"), width=20, heigh=10)
(BNSTMvP | MeAMvP | POAMvP| VMHMvP)/(BNSTMvU | MeAMvU | POAMvU| VMHMv)
dev.off()