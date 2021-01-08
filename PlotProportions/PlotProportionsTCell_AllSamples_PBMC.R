



base = '/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Cluster/PCA40/res3/Plots/'

remove_list = c('HSC','Pro Erythrocyte','Remove','Plasma Cell','GMPC','Pro B Cell','Pre B Cell','CMPC','Erythrocyte','MDPC')
data_run_subset_label = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
data_run_subset_label = data_run_subset_label[,data_run_subset_label$Treatment %in% c('baseline','NBM')]

plot = DimPlot(data_run_subset_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_input  = as.data.frame.matrix(table(data_run_subset_label$sample,Idents(data_run_subset_label) ))


#Import metadata
meta <- data_run_subset_label@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$`Sample Type`[grep("PBMC",meta$`Sample Type`)] <- "PB"
meta$`Sample Type`[grep("Bone Marrow",meta$`Sample Type`)] <- "BM"
meta$Treatment_sampleType = paste0(meta$Treatment,' ', meta$`Sample Type`)
data_input$Group <- meta$Treatment_sampleType[match(rownames(data_input), meta$Sample)]


#Tumor vs Normal comparison
##keep only baseline and NBM samples
##CD14+ to CD16+ Switch
data_input_bl <- data_input[data_input$Group %in% c("baseline BM","baseline PB"),]
data_input_bl$Group <- NULL
##keep only CD14+, CD14+CD16+ and CD16+ data_inputcytes
keep <- c("CD16+ data_input", "SELL+ CD14+ data_input", "CD14+ CD16+ data_input", "MIP1a+ CD14+ data_input", "IFN+ data_input", "TGFb1+ CD14+ data_input")
keep = colnames(data_input_bl)
remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
remove = c('')
data_input_bl_keep <- data_input_bl[,!(colnames(data_input_bl) %in% remove)]
data_input_bl_prop <- data_input_bl_keep/rowSums(data_input_bl_keep)
library(tidyr)
data_input_bl_prop$Patient_name <- rownames(data_input_bl_prop)
data_input_bl_long <- gather(data_input_bl_prop, "Cell_Type","Proportion", -Patient_name)
data_input_bl_long$Type <- "Baseline SMM BM"
data_input_bl_long$Type[grep("PB", data_input_bl_long$Patient_name)] <- "Baseline SMM PB" 

# Get p values for each cell type in all SMM samples
celltype_list = unique(data_input_bl_long$Cell_Type)
for (celltype in celltype_list){
  print(celltype)
  test = wilcox.test(data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) & 
                                               (data_input_bl_long$Type=="Baseline SMM BM")], 
                     data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) &
                                               (data_input_bl_long$Type=="Baseline SMM PB")], conf.int=T)
  print(test$p.value)
  
}

#p-value = 0.0149
tmp1 = data_input_bl_long$Proportion[data_input_bl_long$Type == 'Baseline SMM BM' & 
                                       data_input_bl_long$Cell_Type == 'T Cell']
mean(tmp1)
tmp1 = data_input_bl_long$Proportion[data_input_bl_long$Type == 'Baseline SMM PB' & 
                                       data_input_bl_long$Cell_Type == 'T Cell']
mean(tmp1)
#devtools::install_github("EdwinTh/dutchmasters")
library(dutchmasters)
pdf(paste0(base,"all_prop.pdf"))
ggplot(data_input_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type))  +
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))#+
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
dev.off()

#Volcano plot
volmat <- matrix(nrow=length(unique(data_input_bl_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_SMM", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(data_input_bl_long$Cell_Type)

celltype_list = unique(data_input_bl_long$Cell_Type)
for(ind in 1:length(celltype_list)){
  
  cl <- celltype_list[ind]
  volmat[volmat$Cell_Type == cl,"mean_SMM"] <-   mean(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="Baseline SMM BM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_NBM"] <-   mean(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="Baseline SMM PB"])), na.rm=T)
  
  volmat[volmat$Cell_Type == cl,"Wilcoxon_p"] <- wilcox.test(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & 
                                                                                                                     data_input_bl_long$Type =="Baseline SMM PB"])), 
                                                             as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & 
                                                                                                                     data_input_bl_long$Type =="Baseline SMM BM"])))$p.val
}
volmat$FDR <- p.adjust(as.numeric(as.character(volmat$Wilcoxon_p), method="BH"))
volmat$LFC <- log2(volmat$mean_SMM/volmat$mean_NBM)

library("grid")
crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)

library(ggrepel)
pdf(paste0(base,"Composition Volcano Plot.pdf"))
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin=-3.5, xmax=3.5, ymin=-.2, ymax=3) +
  geom_point() +
  xlim(-3.5,3.5) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")+
  xlab("Log fold-change")+ylab("-log10 p-value") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type))+geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
dev.off()
