# NK
# 1) boxplots for SLAMF7 across my NK cell subtypes (excluding gamma deltas)

gene = 'SLAMF7'
data_input = data_run_subset_label
data_input = data_input[,!(Idents(data_input) %in% c(0:28,'Tgd'))]
Idents(data_input) = as.character(Idents(data_input))

folder = paste0(filepath_cluster,'/Plots/BoxPlots/')
dir.create(folder,recursive = TRUE)
data = as.data.frame( data_input@assays[["RNA"]]@scale.data[gene,])
#data = as.data.frame( data_input@assays[["RNA"]]@data[gene,])

colnames(data) = 'Counts'
data$CellType = as.character(Idents(data_input))

fontSize = 12
pathName = as.character(paste0(folder,'Violin_ ',gene,'_Split_','Cell Type','.pdf') )
pdf(file=pathName,width=14, height=4)
plot = ggplot(data, aes(x=CellType, y=Counts,fill = CellType)) + 
  geom_violin()
plot = plot + ggtitle('SLAMF7 in NK Cell Types') +
  theme_classic()  + theme(
    title =element_text(size=24, color="black"),
    axis.title.x = element_text(color="black", size=fontSize ),
    axis.title.y = element_text(color="black", size=fontSize),
    axis.text= element_text(color="black", size=fontSize),
    legend.text=element_text(size=fontSize),
    legend.title=element_text(size=fontSize),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
print(plot)
dev.off()


# NFkB-high NK cells
# heatmap for NKTR, RUNX3 in NFkB-high cells, 
# comparing baseline to post-treatment (please, collapse C9 & EOT) 
# name the conditions "Baseline" and "Post-treatment".


gene_list = c('NKTR','RUNX3')
data_input = data_run_subset_label
data_input = data_input[,(Idents(data_input) %in% 'NFkB-high')]
Idents(data_input) = as.character(Idents(data_input))
data_input$Treatment[data_input$Treatment == 'baseline'] = 'Baseline'
data_input$Treatment[data_input$Treatment == 'C9D1'] = 'Post-treatment'
data_input$Treatment[data_input$Treatment == 'EOT'] = 'Post-treatment'

data_input = data_input[,data_input$Treatment %in% c('Baseline','Post-treatment')]

folder = paste0(filepath_cluster,'/Plots/Heatmap/')
dir.create(folder,recursive = TRUE)
#data = as.data.frame( data_input@assays[["RNA"]]@scale.data[gene,])

pathName = as.character(paste0(folder,'Heatmap_NK','.pdf') )
pdf(file=pathName,width=8, height=4)
plot = DoHeatmap(data_input, features = gene_list,group.by = 'Treatment')
plot = plot + theme(
  axis.title.x = element_text(color="black", size=12 ),
  axis.title.y = element_text(color="black", size=12),
  axis.text= element_text(color="black", size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12),
  text = element_text(size = 12))
print(plot)
dev.off()


Idents(data_input) = paste0(Idents(data_input),' ',data_input$Treatment)
StackedVlnPlotHelper(data_input,gene_list,
                     folder_heatMap = folder,
                     filename = paste0('Violin_NFkB-highNK.pdf'), width = 10, height = 16)


####################
# cDC2
# 2) heatmap for CD74, HLA-DRB1, HLA-DQA2 for cDC2 comparing baseline to NBM. 
# This is for the paper, so rename baseline to "SMM".


gene_list = c('CD74','HLA-DRB1')
data_input = data_run_subset_label
data_input = data_input[,(Idents(data_input) %in% 'cDC2')]
Idents(data_input) = as.character(Idents(data_input))
data_input$Treatment[data_input$Treatment == 'baseline'] = 'SMM'
data_input = data_input[,data_input$Treatment %in% c('NBM','SMM')]

folder = paste0(filepath_cluster,'/Plots/Heatmap/')
dir.create(folder,recursive = TRUE)
#data = as.data.frame( data_input@assays[["RNA"]]@scale.data[gene,])

pathName = as.character(paste0(folder,'Heatmap_cDC2','.pdf') )
pdf(file=pathName,width=8, height=16)
plot = DoHeatmap(data_input, features = gene_list,group.by = 'Treatment')
plot = plot + theme(
  axis.title.x = element_text(color="black", size=12 ),
  axis.title.y = element_text(color="black", size=12),
  axis.text= element_text(color="black", size=12),
  legend.text=element_text(size=12),
  legend.title=element_text(size=12),
  text = element_text(size = 12))
print(plot)
dev.off()

Idents(data_input) = paste0(Idents(data_input),' ',data_input$Treatment)
StackedVlnPlotHelper(data_input,gene_list,
                     folder_heatMap = folder,
                     filename = paste0('Violin_cDC2.pdf'), width = 10, height = 8*length(gene_list))
