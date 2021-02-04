
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
data_input_raw  = as.data.frame.matrix(table(data_run_subset_label_remove$sample,Idents(data_run_subset_label_remove) ))
#data_input  = as.data.frame.matrix(table(data_run_subset_label_remove$sample,Idents(data_run_subset_label_remove) ))


#Import metadata
meta <- data_run_subset_label_remove@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$Treatment[grep("baseline",meta$Treatment)] <- "Baseline"
data_input = data_input_raw
data_input$Treatment <- meta$Treatment[match(rownames(data_input), meta$Sample)]

data_input$SampleType <- meta$`Sample Type`[match(rownames(data_input), meta$Sample)]

data_input$Dexa <- meta$`Dexa or not`[match(rownames(data_input), meta$Sample)]

data_input$Technology <- meta$Technology[match(rownames(data_input), meta$Sample)]

data_input$Dexa_time = paste0(data_input$Group,' ',data_input$Dexa)

data_input$Treat_SampleType = paste0( data_input$Treatment , ' ', data_input$SampleType)
data_input$Treat_SampleType = gsub('NPBMC NPBMC','NPBMC',data_input$Treat_SampleType )


data_input$Tech_Treat_SampleType = paste0(data_input$Technology,' ', data_input$Treatment , ' ', data_input$SampleType)
data_input$Tech_Treat_SampleType = gsub('NPBMC NPBMC','NPBMC',data_input$Tech_Treat_SampleType )

data_input$Group = data_input$Tech_Treat_SampleType 

#####################################################
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

celltype = 'Mono'
celltype = 'DC'

#Volcano plot NPBMC Vs PBMC
if (celltype == 'Mono'){
  
  Ident_order = c('SELL+ CD14+ Mono','TGFb1+ CD14+ Mono',
                'IFN+ Mono','CD14+ CD16+ Mono','CD16+ Mono')
}
if (celltype == 'DC'){
  Ident_order = c('cDC1','cDC2','prDC')
  
}
#remove = c('0','11','12','17','18','21','Erythrocyte','T/NK Doublet')
data_input_keep = as.matrix(data_input[Ident_order])
data_input_keep = as.data.frame(data_input_keep)
#colnames(data_input_keep) = Ident_order
prop <- data_input_keep/rowSums(data_input_keep) # Compare to all cell types


prop_meta = prop
prop_meta$Patient_name <- rownames(prop)

##
prop_meta$Group <- data_input$Group
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)
prop_meta_long$Cell_Type = factor(prop_meta_long$Cell_Type, levels = Ident_order)

unique(prop_meta_long$Group)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Microwell NPBMC','10x NPBMC'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'Microwell NPBMC',var2 = '10x NPBMC',
                  base = paste0(base,celltype,'/'), xmin = -3,xmax = 3)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('10x Baseline PBMC','10x NPBMC'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = '10x Baseline PBMC',var2 = '10x NPBMC',
                  base = paste0(base,celltype,'/'), xmin = -3,xmax = 3)


prop_meta$Group <- data_input$Treat_SampleType
prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group)
prop_meta_long$Cell_Type = factor(prop_meta_long$Cell_Type, levels = Ident_order)

unique(prop_meta_long$Group)
data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline Bone Marrow','Baseline PBMC'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'Baseline Bone Marrow',var2 = 'Baseline PBMC',
                  base = paste0(base,celltype,'/'), xmin = -3,xmax = 3)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Group %in% c('Baseline Bone Marrow','NPBMC'),]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'Baseline Bone Marrow',var2 = 'NPBMC',
                  base = paste0(base,celltype,'/'), xmin = -3,xmax = 3)

