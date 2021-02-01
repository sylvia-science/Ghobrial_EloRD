
library(tidyr)
library(dutchmasters)
library(ggrepel)
library("grid")

source('~/Desktop/scRNA/Code/Plot_func.R')


base = paste0(filepath_cluster,'Plots/')

#base = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/',
#              'Plots/Paper/Proportions/')
dir.create(base)

Idents(data_harmony_run_label_remove)[Idents(data_harmony_run_label_remove) == 'CD8+ T-cell'] = 'T-cell'
data_input_raw  = as.data.frame.matrix(table(data_harmony_run_label$sample,Idents(data_harmony_run_label) ))
#data_input  = as.data.frame.matrix(table(data_run_subset_label_remove$sample,Idents(data_run_subset_label_remove) ))


#Import metadata
meta <- data_harmony_run_label@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$Treatment[grep("baseline",meta$Treatment)] <- "Baseline"
data_input = data_input_raw
data_input$Group <- meta$Treatment[match(rownames(data_input), meta$Sample)]

data_input$Dexa <- meta$`Dexa or not`[match(rownames(data_input), meta$Sample)]

data_input$SampleType <- meta$`Sample Type`[match(rownames(data_input), meta$Sample)]

data_input$Dexa_time = paste0(data_input$Group,' ',data_input$Dexa)

#Tumor vs Normal comparison
##keep only baseline and NBM samples
##CD14+ to CD16+ Switch
data_input_bl <- data_input[data_input$Group %in% c("Baseline","NBM"),]
data_input_bl$Group <- NULL
data_input_bl$Dexa <- NULL
data_input_bl$Dexa_time = NULL

write.csv(data_input_bl,paste0(base,'cellCountsAll.csv'))
##keep only CD14+, CD14+CD16+ and CD16+ data_inputcytes
keep <- c("CD16+ data_input", "SELL+ CD14+ data_input", "CD14+ CD16+ data_input", "MIP1a+ CD14+ data_input", "IFN+ data_input", "TGFb1+ CD14+ data_input")

keep = c('cDC1','cDC2','sDC','prDC','dDC')
#data_input_bl_keep <- data_input_bl[,(colnames(data_input_bl) %in% keep)]

remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
data_input_bl_keep <- data_input_bl[,!(colnames(data_input_bl) %in% remove)]

data_input_bl_prop <- data_input_bl_keep/rowSums(data_input_bl_keep)


data_input_bl_prop$Patient_name <- rownames(data_input_bl_prop)
data_input_bl_long <- gather(data_input_bl_prop, "Cell_Type","Proportion", -Patient_name)
data_input_bl_long$Type <- "SMM"
data_input_bl_long$Type[grep("NBM", data_input_bl_long$Patient_name)] <- "NBM" 

tmp = data_input_bl_long[data_input_bl_long$Cell_Type == 'T-cell',]
# Get p values for each cell type in all SMM samples
celltype_list = unique(data_input_bl_long$Cell_Type)
# for (celltype in celltype_list){
#   print(celltype)
#   var1_prop = data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) & 
#                                               (data_input_bl_long$Type=="NBM")]
#   var2_prop = data_input_bl_long$Proportion[(data_input_bl_long$Cell_Type==celltype) & 
#                                               (data_input_bl_long$Type=="SMM")]
#   test = wilcox.test(var1_prop, var2_prop, conf.int=T)
#   print(test$p.value)
#   
# }


pdf(paste0(base,celltype,".pdf"))
plot = ggplot(data_input_bl_long) + geom_boxplot(aes(x=Type,y=Proportion, fill=Cell_Type))  +
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))#+
#theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
print(plot)
dev.off()





#Volcano plot Baseline Vs EOT
celltype = 'Mono'
celltype = 'DC'

if (celltype == 'Mono'){
  Ident_order = c('SELL+ CD14+ Mono','sMono','MIP1a+ CD14+ Mono','TGFb1+ CD14+ Mono',
                  'sCD14+ Mono','IFN+ Mono','CD14+ CD16+ Mono','CD16+ Mono')
  
}
if (celltype == 'DC'){
  Ident_order = c('cDC1','cDC2','prDC', 'sDC')
  
}

Ident_order = c('Pro-B-cell','Pre B-cell','B-cell')

#remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
data_input_bl_keep <- data_input[,colnames(data_input) %in% Ident_order]
prop <- data_input_bl_keep/rowSums(data_input_raw) # Compare to all cell types

data_input_bl_keep$Group <- NULL
data_input_bl_keep$Dexa <- NULL
data_input_bl_keep$Dexa_time <- NULL
#prop <- data_input_bl_keep/rowSums(data_input_bl_keep)


prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input$Group
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)
prop_meta_long$Cell_Type = factor(prop_meta_long$Cell_Type, levels = Ident_order)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline','EOT'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'EOT',var2 = 'Baseline',base)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline','C9D1'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'C9D1',var2  = 'Baseline',base)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('NBM','Baseline'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'Baseline',var2  = 'NBM',base, str = paste0('_',celltype))

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('NBM','EOT'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1 = 'EOT',var2  = 'NBM',base)


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


##############


Ident_order = c('B-cell')

#remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
data_input_keep = as.matrix(data_input[Ident_order])
data_input_keep = as.data.frame(data_input_keep)
colnames(data_input_keep) = Ident_order
prop <- data_input_keep/rowSums(data_input_raw) # Compare to all cell types


prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
Treatment = meta$Treatment[match(rownames(data_input), meta$Sample)]
sample_type = meta$`Sample Type`[match(rownames(data_input), meta$Sample)]

prop_meta$Group <- data_input$SampleType
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)
prop_meta_long$Cell_Type = factor(prop_meta_long$Cell_Type, levels = Ident_order)


data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('NPBMC','PBMC'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'NPBMC',var2 = 'PBMC',base)
