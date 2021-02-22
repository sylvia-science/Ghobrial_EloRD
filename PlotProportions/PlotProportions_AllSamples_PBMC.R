


source('~/Desktop/scRNA/Code/Plot_func.R')
base = '/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Cluster/PCA40/res3/Plots/'

remove_list = c('HSC','Pro Erythrocyte','Plasma Cell','GMPC',
                'Pro B-cell','Pre B-cell','CMPC','Erythrocyte','MDPC', 1:20)
data_run_label = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
data_run_label = data_harmony_run_label
data_run_label = data_run_label[,data_run_label$Treatment %in% c('baseline','NBM')]

plot = DimPlot(data_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_input  = as.data.frame.matrix(table(data_run_label$sample,Idents(data_run_label) ))


str = '_keepAll'
#Import metadata
meta <- data_run_label@meta.data
meta$Treatment[grep("EOT",meta$Treatment)] <- "EOT"
meta$`Sample Type`[grep("PBMC",meta$`Sample Type`)] <- "PB"
meta$`Sample Type`[grep("Bone Marrow",meta$`Sample Type`)] <- "BM"
meta$Treatment[grep("baseline",meta$Treatment)] <- "Baseline"

meta$Treatment_sampleType = paste0(meta$Treatment,' ', meta$`Sample Type`)

data_input$Treatment <- meta$Treatment[match(rownames(data_input), meta$Sample)]

data_input$Group <- meta$Treatment[match(rownames(data_input), meta$Sample)]

data_input$Dexa <- meta$`Dexa or not`[match(rownames(data_input), meta$Sample)]

data_input$Dexa_time = paste0(data_input$Group,' ',data_input$Dexa)

data_input$sample_type = meta$'Sample Type'[match(rownames(data_input), meta$Sample)]

remove = c('')
data_input_bl_keep <- data_input[,!(colnames(data_input) %in% remove)]
data_input_bl_keep$Treatment <- NULL
data_input_bl_keep$Group <- NULL
data_input_bl_keep$Dexa <- NULL
data_input_bl_keep$Dexa_time <- NULL
data_input_bl_keep$sample_type <- NULL
prop <- data_input_bl_keep/rowSums(data_input_bl_keep)

prop_meta = prop
prop_meta$Patient_name <- rownames(prop)
prop_meta$Group <- data_input$sample_type
prop_meta$Treatment <- data_input$Treatment

prop_meta_long <- gather(prop_meta, "Cell_Type","Proportion", -Patient_name,-Group,-Treatment)

data_input_bl_long_subset = prop_meta_long[prop_meta_long$Treatment == 'Baseline',]
VolcanoPlotHelper(data_input_bl_long_subset,var1  = 'BM',var2 = 'PB',paste0(base,str))


####

#devtools::install_github("EdwinTh/dutchmasters")
library(dutchmasters)
pdf(paste0(base,str,"all_prop_BMVsPB.pdf"))
ggplot(data_input_bl_long_subset) + geom_boxplot(aes(x=Group,y=Proportion, fill=Cell_Type))  +
  #scale_fill_manual(values=c("#D5BF98","#AF7366","#8B6C4F","#CDD4E4","#E3C78F","#78A8D1"))+ 
  theme_bw() + xlab("") + theme(axis.text.x=element_text(size=12))#+
  #theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.box.margin=margin(c(1,1,1,1)))
dev.off()

