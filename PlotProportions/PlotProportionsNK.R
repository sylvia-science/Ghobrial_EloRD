

  
base = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/NK/Cluster/PCA30/res3/Plots/'
mono  = as.data.frame.matrix(table(data_run_subset_label$sample,Idents(data_run_subset_label) ))

#Import metadata
meta <- data_run_subset_label@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
mono$Group <- meta$Treatment[match(rownames(mono), meta$Sample)]


#Tumor vs Normal comparison
##keep only baseline and NBM samples
##CD14+ to CD16+ Switch
mono_bl <- mono[mono$Group %in% c("baseline","NBM"),]
mono_bl$Group <- NULL
##keep only CD14+, CD14+CD16+ and CD16+ Monocytes
keep <- c("CD16+ Mono", "SELL+ CD14+ Mono", "CD14+ CD16+ Mono", "MIP1a+ CD14+ Mono", "IFN+ Mono", "TGFb1+ CD14+ Mono")
keep = colnames(mono_bl)
remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
mono_bl_keep <- mono_bl[,!(colnames(mono_bl) %in% remove)]
mono_bl_prop <- mono_bl_keep/rowSums(mono_bl_keep)
library(tidyr)
mono_bl_prop$Patient_name <- rownames(mono_bl_prop)
mono_bl_long <- gather(mono_bl_prop, "Cell_Type","Proportion", -Patient_name)
mono_bl_long$Type <- "SMM"
mono_bl_long$Type[grep("NBM", mono_bl_long$Patient_name)] <- "NBM" 

# Get p values for each cell type in all SMM samples
celltype_list = unique(mono_bl_long$Cell_Type)
for (celltype in celltype_list){
  print(celltype)
  test = wilcox.test(mono_bl_long$Proportion[(mono_bl_long$Cell_Type==celltype) & 
                                               (mono_bl_long$Type=="NBM")], 
                     mono_bl_long$Proportion[(mono_bl_long$Cell_Type==celltype) &
                                               (mono_bl_long$Type=="SMM")], conf.int=T)
  print(test$p.value)
  
}

#p-value = 0.0149

#devtools::install_github("EdwinTh/dutchmasters")
library(dutchmasters)
pdf(paste0(base,"all_prop.pdf"))
ggplot(mono_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type))  +
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
dev.off()

#Volcano plot
volmat <- matrix(nrow=length(unique(mono_bl_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_SMM", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(mono_bl_long$Cell_Type)
for(ind in 1:length(unique(mono_bl_long$Cell_Type))){
  cl <- unique(mono_bl_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl,"mean_SMM"] <-   mean(as.numeric(as.character(mono_bl_long$Proportion[mono_bl_long$Cell_Type==cl & mono_bl_long$Type =="SMM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_NBM"] <-   mean(as.numeric(as.character(mono_bl_long$Proportion[mono_bl_long$Cell_Type==cl & mono_bl_long$Type =="NBM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"Wilcoxon_p"] <- wilcox.test(as.numeric(as.character(mono_bl_long$Proportion[mono_bl_long$Cell_Type==cl & mono_bl_long$Type =="NBM"])), as.numeric(as.character(mono_bl_long$Proportion[mono_bl_long$Cell_Type==cl & mono_bl_long$Type =="SMM"])))$p.val
}
volmat$FDR <- p.adjust(as.numeric(as.character(volmat$Wilcoxon_p), method="BH"))
volmat$LFC <- log2(volmat$mean_SMM/volmat$mean_NBM)

library("grid")
crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)

library(ggrepel)
pdf(paste0(base,"Composition Volcano Plot.pdf"))
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin=-3.5, xmax=3.5, ymin=-.2, ymax=2.5) +
  geom_point() +
  xlim(-3.5,3.5) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")+
  xlab("Log fold-change")+ylab("-log10 p-value") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type))+geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
dev.off()


keep <- c("SELL+ CD14+ Mono", "TGFb1+ CD14+ Mono")
mono_bl_keep <- mono_bl[,colnames(mono_bl) %in% keep]
mono_bl_prop <- mono_bl_keep/rowSums(mono_bl_keep)
library(tidyr)
mono_bl_prop$Patient_name <- rownames(mono_bl_prop)
mono_bl_long <- gather(mono_bl_prop, "Cell_Type","Proportion", -Patient_name)
mono_bl_long$Type <- "SMM"
mono_bl_long$Type[grep("NBM", mono_bl_long$Patient_name)] <- "NBM" 

wilcox.test(mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="SELL+ CD14+ Mono") & (mono_bl_long$Type=="SMM")], mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="TGFb1+ CD14+ Mono") & (mono_bl_long$Type=="SMM")], conf.int=T)
#p-value = 0.909

wilcox.test(mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="SELL+ CD14+ Mono") & (mono_bl_long$Type=="NBM")], mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="TGFb1+ CD14+ Mono") & (mono_bl_long$Type=="NBM")], conf.int=T)
#p-value = 0.02101

wilcox.test(mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="SELL+ CD14+ Mono") & (mono_bl_long$Type=="NBM")], mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="SELL+ CD14+ Mono") & (mono_bl_long$Type=="SMM")], conf.int=T)
#p-value = 0.2697

wilcox.test(mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="TGFb1+ CD14+ Mono") & (mono_bl_long$Type=="NBM")], mono_bl_long$Proportion[(mono_bl_long$Cell_Type=="TGFb1+ CD14+ Mono") & (mono_bl_long$Type=="SMM")], conf.int=T)
#p-value = 0.2697

pdf(paste0(base,"Program Switch.pdf"))
ggplot(mono_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type)) + annotate("text", x= 1, y=0.1, label="p-value = 0.02", size = 5) + annotate("text", x= 2, y = 1.05, label = "p-value = 0.9", size = 5) + scale_fill_manual(values=c("#E3C78F","#78A8D1"))+ theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))
dev.off()


###Monocyte proportions over time with treatment
mono_bl <- mono
mono_bl$Group <- NULL
##keep only CD14+, CD14+CD16+ and CD16+ Monocytes
keep <- c("CD16+ Mono", "SELL+ CD14+ Mono", "CD14+ CD16+ Mono", "MIP1a+ CD14+ Mono", "IFN+ Mono", "TGFb1+ CD14+ Mono")
remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')

mono_bl_keep <- mono_bl[,!(colnames(mono_bl) %in% remove)]
mono_bl_prop <- mono_bl_keep/rowSums(mono_bl_keep)
library(tidyr)
mono_bl_prop$Patient_name <- rownames(mono_bl_prop)
mono_bl_long <- gather(mono_bl_prop, "Cell_Type","Proportion", -Patient_name)
mono_bl_long$Type <- factor(meta$Treatment[match(mono_bl_long$Patient_name, meta$Sample)], levels=c("NBM", "baseline","C9D1", "EOT"))

pdf(paste0(base,"Treatment effect on NK Proportions.pdf"))
ggplot(mono_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type))  + 
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + 
  theme(axis.text.x=element_text(size=12))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
dev.off()

####Plot boxplots of proportions with lines connecting individual patients' dots
mono_bl_long$Patient <- meta$`Patient Number`[match(mono_bl_long$Patient_name, meta$Sample)]
mono_bl_long$BOR <- meta$Best_Overall_Response[match(mono_bl_long$Patient_name, meta$Sample)]

library(wesanderson)
pdf(paste0(base,"Lineplots Per NK Subtype.pdf"))
for(n in 1:length(unique(mono_bl_long$Cell_Type))){
  cell_type <- as.character(unique(mono_bl_long$Cell_Type))[n]
  line_input <- mono_bl_long[(mono_bl_long$Cell_Type ==cell_type) & (!mono_bl_long$Type %in% "NBM"),]
  print(ggplot(line_input, aes(x = Type, y = Proportion)) +
          geom_boxplot() +
          geom_point(color="black", size=2) +
          geom_line(aes(group=Patient, color=BOR), size = 1) +
          theme_classic() + scale_color_manual(values=wes_palette("FantasticFox1"), 
                                               labels=c("sCR","CR","VGPR","PR","MR")) + xlab("") +
          theme(axis.text.x=element_text(size=12)) + ylab(paste0(cell_type," Proportion")))
}
dev.off()


######Survival analysis
#Import survival data & create survival object
library(survival)
surv <- read.csv("Data/Data_From_Rob/14338-manu.csv")
surv$SurvObj <- with(surv, Surv(pfs_time, pfs_ind == 1))
mono_bl <- mono[mono$Group %in% c("baseline"),]
mono_bl$Group <- NULL
##keep only CD14+, CD14+CD16+ and CD16+ Monocytes
keep <- c("CD16+ Mono", "SELL+ CD14+ Mono", "CD14+ CD16+ Mono", "MIP1a+ CD14+ Mono", "IFN+ Mono", "TGFb1+ CD14+ Mono")
mono_bl_keep <- mono_bl[,colnames(mono_bl) %in% keep]
mono_bl_prop <- mono_bl_keep/rowSums(mono_bl_keep)
mono_bl_prop$casenum <- meta$Patient.Number[match(rownames(mono_bl_prop),meta$Sample)]
mono_surv <- merge(mono_bl_prop, surv, by="casenum")
colnames(mono_surv)[2:7] <- c("CD16_Mono","TGFb1_CD14_Mono","SELL_CD14_Mono","CD14_CD16_Mono","IFN_Mono","MIP1a_CD14_Mono")
#Fit the Cox model & plot
fit <- coxph(SurvObj ~ CD16_Mono + TGFb1_CD14_Mono + SELL_CD14_Mono + CD14_CD16_Mono + IFN_Mono + MIP1a_CD14_Mono,  data = mono_surv)
summary(fit)

#Retry with CD14, CD14/CD16 and CD16 only (based on a result showing S1 is enriched in progressors)
mono_bl_keep <- mono_bl[,colnames(mono_bl) %in% keep]
#mono_bl_keep$`CD14+ Mono` <- mono_bl_keep$`TGFb1+ CD14+ Mono` + mono_bl_keep$`SELL+ CD14+ Mono` + mono_bl_keep$`IFN+ Mono` + mono_bl_keep$`MIP1a+ CD14+ Mono`
#mono_bl_keep <- mono_bl_keep[,c("CD14+ Mono","CD14+ CD16+ Mono","CD16+ Mono")]
mono_bl_prop <- mono_bl_keep/rowSums(mono_bl_keep)
mono_bl_prop$casenum <- meta$Patient.Number[match(rownames(mono_bl_prop),meta$Sample)]
mono_surv <- merge(mono_bl_prop, surv, by="casenum")
#colnames(mono_surv)[2:4] <- c("CD14_Mono","CD14_CD16_Mono","CD16_Mono")
colnames(mono_surv)[2:7] <- c("CD16_Mono", "TGFb1_CD14_Mono","SELL_CD14_Mono","CD14_CD16_Mono","IFN_Mono","MIP1a_CD14_Mono")
#mono_surv <- mono_surv[,c("CD14_Mono","CD14_CD16_Mono","CD16_Mono","pfs_ind", "pfs_time")]
mono_surv <- mono_surv[,c("CD16_Mono", "TGFb1_CD14_Mono","SELL_CD14_Mono","CD14_CD16_Mono","IFN_Mono","MIP1a_CD14_Mono", "casenum")]

mono_surv_long <- gather(mono_surv,"Cell_Type","Proportion",-casenum)
ggplot(mono_surv_long,aes(factor(casenum),Proportion, fill= Cell_Type)) + geom_bar(position="fill", stat="identity")


pdf(paste0(base,"Lineplots Per Monocyte Subtype by PFS.pdf"))
for(n in 1:length(unique(mono_bl_long$Cell_Type))){
  cell_type <- as.character(unique(mono_bl_long$Cell_Type))[n]
  line_input <- mono_bl_long[(mono_bl_long$Cell_Type ==cell_type) & (!mono_bl_long$Type %in% "NBM"),]
  line_input$pfs_time <- surv$pfs_time[match(line_input$Patient, surv$casenum)]
  print(ggplot(line_input, aes(x = Type, y = Proportion)) +
          geom_boxplot() +
          geom_point(color="black", size=2) +
          geom_line(aes(group=Patient, color=pfs_time), size = 1) +
          scale_color_gradient(low="red3",high="cornsilk2")+
          theme_classic() + xlab("") +theme(axis.text.x=element_text(size=12)) + ylab(paste0(cell_type," Proportion")))
}
dev.off()

  
  
  
  