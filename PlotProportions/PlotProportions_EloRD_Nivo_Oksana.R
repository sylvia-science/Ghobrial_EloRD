

library(tidyr)
library(dutchmasters)

source('~/Desktop/scRNA/Code/Plot_func.R')

celltype_orig = celltype
base = '/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Cluster/PCA40/res3/Plots/'
base = paste0(filepath_cluster,'/Plots/')
dir.create(base)
unique(Idents(data_run_subset_label))[18]

remove_list = c('HSC','Pro Erythrocyte','Plasma Cell','GMPC',
                'Pro B-cell','Pre B-cell','CMPC','Erythrocyte','MDPC','Remove', 'dMono','dNK')

data_run_label = data_run_subset_label[,!(Idents(data_run_subset_label) %in% remove_list)]
#data_run_label = data_run_subset_label
unique(Idents(data_run_label))
#data_run_label = data_run_label[,data_run_label$Treatment %in% c('Pre-treatment','NBM')]

plot = DimPlot(data_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_input  = as.data.frame.matrix(table(data_run_label$sample,Idents(data_run_label) ))


str = '_keepAll'
meta <- data_run_label@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$`Sample Type`[grep("Bone Marrow",meta$`Sample Type`)] <- "BM"
meta$Treatment[grep("baseline",meta$Treatment)] <- "Baseline"

meta$Treatment_sampleType = paste0(meta$Treatment,' ', meta$`Sample Type`)

data_input$Treatment <- meta$Treatment[match(rownames(data_input), meta$Sample)]
data_input$Group <- meta$Treatment[match(rownames(data_input), meta$Sample)]
data_input$Dexa <- meta$`Dexa or not`[match(rownames(data_input), meta$Sample)]
data_input$Dexa_time = paste0(data_input$Group,' ',data_input$Dexa)
data_input$sample_type = meta$'Sample Type'[match(rownames(data_input), meta$Sample)]
data_input$Diagnosis =  meta$Diagnosis[match(rownames(data_input), meta$Sample)]
data_input$Response =  meta$Response[match(rownames(data_input), meta$Sample)]
data_input$Study =  meta$Study[match(rownames(data_input), meta$Sample)]

#########################################
## Plot Diagnosis in pre treatment
#########################################
data_input_subset = data_input[data_input$Treatment %in% c('Pre-treatment',"Pre-treatment (collapse with old sample)",'NBM'),]


remove = c('')
data_input_keep <- data_input_subset[,!(colnames(data_input_subset) %in% remove)]

data_input_keep$Treatment <- NULL
data_input_keep$Group <- NULL
data_input_keep$Dexa <- NULL
data_input_keep$Dexa_time <- NULL
data_input_keep$sample_type <- NULL
data_input_keep$Diagnosis = NULL
data_input_keep$Response = NULL
data_input_keep$Study = NULL

prop <- data_input_keep/rowSums(data_input_keep)

prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input_subset$sample_type
prop_meta$Treatment <- data_input_subset$Treatment
prop_meta$Diagnosis <- data_input_subset$Diagnosis

prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group,-Treatment,-Diagnosis)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Treatment == 'Baseline',]
#VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'BM',var2 = 'PB',paste0(base,str))

prop_meta_long$Diagnosis = factor(prop_meta_long$Diagnosis, levels = c('NBM','MGUS','Low Risk SMM','High Risk SMM','Newly Diagnosed MM'))

####
library(dutchmasters)

str = ''
celltype_list = unique(prop_meta_long$Cell_Type)
celltype = celltype_list[1]
dir.create(paste0(base,'Boxplot/'))
pdf(paste0(base,'Boxplot/',celltype_orig,"_Diagnosis",".pdf"),width = 12)

for (celltype in celltype_list){

  prop_meta_long_celltype = prop_meta_long[prop_meta_long$Cell_Type == celltype,]
  plot = ggplot(prop_meta_long_celltype) + geom_boxplot(aes(x=Diagnosis,y=Proportion, fill=Diagnosis))  +
    #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
    theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))+ 
    ggtitle(paste0(celltype)) + theme(plot.title=element_text(hjust=0.5))
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
  print(plot)
  

}

dev.off()


#############################################
## Plot response in pre treatment in Nivo
#############################################
data_input_subset = data_input[data_input$Treatment %in% c('Pre-treatment','NBM'),]

data_input_subset = data_input_subset[data_input_subset$Study %in% c('Nivolumab Trial','NBM'),]
data_input_subset = data_input_subset[data_input_subset$Response != '',]
remove = c('')
data_input_keep <- data_input_subset[,!(colnames(data_input_subset) %in% remove)]

data_input_keep$Treatment <- NULL
data_input_keep$Group <- NULL
data_input_keep$Dexa <- NULL
data_input_keep$Dexa_time <- NULL
data_input_keep$sample_type <- NULL
data_input_keep$Diagnosis = NULL
data_input_keep$Response = NULL
data_input_keep$Study = NULL

prop <- data_input_keep/rowSums(data_input_keep)

prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input_subset$sample_type
prop_meta$Treatment <- data_input_subset$Treatment
prop_meta$Diagnosis <- data_input_subset$Diagnosis
prop_meta$Response <- data_input_subset$Response

prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group,-Treatment,-Diagnosis,-Response)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Treatment == 'Baseline',]
#VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'BM',var2 = 'PB',paste0(base,str))

prop_meta_long$Response = factor(prop_meta_long$Response, levels = c('Healthy','Stable Disease','PR then Progression','PR','VGPR','CR'))

####

str = ''
celltype_list = unique(prop_meta_long$Cell_Type)
celltype = celltype_list[1]
dir.create(paste0(base,'Boxplot/'))
pdf(paste0(base,'Boxplot/',celltype_orig, str,"_PreTreatmentResponse",".pdf"),width = 12)

for (celltype in celltype_list){
  
  prop_meta_long_celltype = prop_meta_long[prop_meta_long$Cell_Type == celltype,]
  plot = ggplot(prop_meta_long_celltype) + geom_boxplot(aes(x=Response,y=Proportion, fill=Response))  +
    #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
    theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12)) + 
  ggtitle(paste0(celltype)) + theme(plot.title=element_text(hjust=0.5))
  
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
  print(plot)
  
  
}

dev.off()

##



#############################################
## Plot Pre/post treatment in Nivo trial
#############################################
data_input_subset = data_input[data_input$Treatment %in% c('Pre-treatment','Post-treatment','NBM'),]

data_input_subset = data_input_subset[data_input_subset$Study %in% c('Nivolumab Trial','NBM'),]
#data_input_subset = data_input_subset[data_input_subset$Response != '',]
remove = c('')
data_input_keep <- data_input_subset[,!(colnames(data_input_subset) %in% remove)]

data_input_keep$Treatment <- NULL
data_input_keep$Group <- NULL
data_input_keep$Dexa <- NULL
data_input_keep$Dexa_time <- NULL
data_input_keep$sample_type <- NULL
data_input_keep$Diagnosis = NULL
data_input_keep$Response = NULL
data_input_keep$Study = NULL

prop <- data_input_keep/rowSums(data_input_keep)

prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input_subset$sample_type
prop_meta$Treatment <- data_input_subset$Treatment
prop_meta$Diagnosis <- data_input_subset$Diagnosis
prop_meta$Response <- data_input_subset$Response

prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group,-Treatment,-Diagnosis,-Response)


prop_meta_long$Response = factor(prop_meta_long$Response, levels = c('Healthy','Stable Disease','PR then Progression','PR','VGPR','CR'))
prop_meta_long$Treatment = factor(prop_meta_long$Treatment, levels = c('NBM','Pre-treatment','Post-treatment'))


####

str = ''
celltype_list = unique(prop_meta_long$Cell_Type)
celltype = celltype_list[1]
dir.create(paste0(base,'Boxplot/'))
pdf(paste0(base,'Boxplot/',celltype_orig, str,"_NivoTreatment",".pdf"),width = 12)

for (celltype in celltype_list){
  
  prop_meta_long_celltype = prop_meta_long[prop_meta_long$Cell_Type == celltype,]
  plot = ggplot(prop_meta_long_celltype) + geom_boxplot(aes(x=Treatment,y=Proportion, fill=Treatment))  +
    #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
    theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12)) + 
    ggtitle(paste0(celltype)) + theme(plot.title=element_text(hjust=0.5))
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
  print(plot)
  
}

dev.off()

data_input_volcano_long_subset = prop_meta_long[prop_meta_long$Treatment %in% c('Pre-treatment','Post-treatment'),]
data_input_volcano_long_subset$Group = data_input_volcano_long_subset$Treatment
VolcanoPlotHelper(data_input_volcano_long_subset,var1  = 'Post-treatment',var2 = 'Pre-treatment',paste0(base,'Nivo_',celltype_orig,'_'))

##



