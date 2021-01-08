

source('~/Desktop/scRNA/Code/Plot_func.R')
base = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/T Cell/Cluster/PCA30/res3.5/Plots/'
base = paste0(filepath_cluster,'Plots/')
dir.create(base)
data_input  = as.data.frame.matrix(table(data_run_subset_label_remove$sample,Idents(data_run_subset_label_remove) ))

#Import metadata
meta <- data_run_subset_label_remove@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$Treatment[grep("baseline",meta$Treatment)] <- "Baseline"

data_input$Group <- meta$Treatment[match(rownames(data_input), meta$Sample)]

data_input$Dexa <- meta$`Dexa or not`[match(rownames(data_input), meta$Sample)]

data_input$Dexa_time = paste0(data_input$Group,' ',data_input$Dexa)

#Tumor vs Normal comparison
##keep only baseline and NBM samples
##CD14+ to CD16+ Switch
data_input_bl <- data_input[data_input$Group %in% c("Baseline","NBM"),]
data_input_bl$Group <- NULL
data_input_bl$Dexa <- NULL
data_input_bl$Dexa_time <- NULL
##keep only CD14+, CD14+CD16+ and CD16+ data_inputcytes
keep <- c("CD16+ data_input", "SELL+ CD14+ data_input", "CD14+ CD16+ data_input", "MIP1a+ CD14+ data_input", "IFN+ data_input", "TGFb1+ CD14+ data_input")
keep = colnames(data_input_bl)
remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
data_input_bl_keep <- data_input_bl[,!(colnames(data_input_bl) %in% remove)]
data_input_bl_prop <- data_input_bl_keep/rowSums(data_input_bl_keep)
library(tidyr)
data_input_bl_prop$Patient_name <- rownames(data_input_bl_prop)
data_input_bl_long <- gather(data_input_bl_prop, "Cell_Type","Proportion", -Patient_name)
data_input_bl_long$Type <- "SMM"
data_input_bl_long$Type[grep("NBM", data_input_bl_long$Patient_name)] <- "NBM" 

# Get p values for each cell type in all SMM samples
celltype_list = unique(data_input_bl_long$Cell_Type)
for (celltype in celltype_list){
  print(celltype)
  var1_prop = data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) & 
                                              (data_input_bl_long$Type=="NBM")]
  var2_prop = data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) & 
                                              (data_input_bl_long$Type=="SMM")]
  test = wilcox.test(var1_prop, var2_prop, conf.int=T)
  print(test$p.value)
  
}

#p-value = 0.0149

#devtools::install_github("EdwinTh/dutchmasters")
library(dutchmasters)
pdf(paste0(base,"all_prop.pdf"))
ggplot(data_input_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type))  +
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))+
  theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
dev.off()

#Volcano plot
volmat <- matrix(nrow=length(unique(data_input_bl_long$Cell_Type)), ncol=6)
colnames(volmat) <- c("Cell_Type","Wilcoxon_p", "FDR", "mean_SMM", "mean_NBM", "LFC")
volmat <- data.frame(volmat)
volmat$Cell_Type <- unique(data_input_bl_long$Cell_Type)
for(ind in 1:length(unique(data_input_bl_long$Cell_Type))){
  cl <- unique(data_input_bl_long$Cell_Type)[ind]
  volmat[volmat$Cell_Type == cl,"mean_SMM"] <-   mean(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="SMM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"mean_NBM"] <-   mean(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="NBM"])), na.rm=T)
  volmat[volmat$Cell_Type == cl,"Wilcoxon_p"] <- wilcox.test(as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="NBM"])), as.numeric(as.character(data_input_bl_long$Proportion[data_input_bl_long$Cell_Type==cl & data_input_bl_long$Type =="SMM"])))$p.val
}
volmat$FDR <- p.adjust(as.numeric(as.character(volmat$Wilcoxon_p), method="BH"))
volmat$LFC <- log2(volmat$mean_SMM/volmat$mean_NBM)

write.csv(volmat,paste0(base,"volcano_data_NBMVsSMM.csv"))


library("grid")
crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)

library(ggrepel)
pdf(paste0(base,"Composition Volcano Plot_NBMVsSMM.pdf"))
ggplot(volmat,aes(x=LFC, y=-log10(Wilcoxon_p)), size=4) +
  annotation_custom(g, xmin=-3.5, xmax=3.5, ymin=-.2, ymax=3) +
  geom_point() +
  xlim(-3.5,3.5) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")+
  xlab("Log fold-change")+ylab("-log10 p-value") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(aes(x=LFC, y=-log10(Wilcoxon_p), label=Cell_Type))+geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
dev.off()






#Volcano plot Baseline Vs EOT

remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
remove = c('')
data_input_bl_keep <- data_input[,!(colnames(data_input) %in% remove)]
data_input_bl_keep$Group <- NULL
data_input_bl_keep$Dexa <- NULL
data_input_bl_keep$Dexa_time <- NULL
prop <- data_input_bl_keep/rowSums(data_input_bl_keep)

prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input$Group
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline','EOT'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'EOT',var2 = 'Baseline')

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline','C9D1'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'C9D1',var2  = 'Baseline')

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('NBM','Baseline'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'Baseline',var2  = 'NBM')

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('NBM','EOT'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'EOT',var2  = 'NBM')


prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <-data_input$Dexa_time
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('C9D1 Yes','C9D1 No'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'C9D1 Yes',var2  = 'C9D1 No')

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('EOT Yes','EOT No'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'EOT Yes',var2  = 'EOT No')

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline Yes','Baseline No'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'Baseline Yes',var2  = 'Baseline No')
