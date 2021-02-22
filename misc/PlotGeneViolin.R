data_input = data_run_subset_label_remove

gene = 'B3GAT1'
data_subset_gene = t(data_input[gene,][["RNA"]]@data)
data_subset_gene = as.data.frame(data_subset_gene)
colnames(data_subset_gene) = 'gene'
data_subset_gene$label = Idents(data_input)
data_subset_gene$sample = data_input$sample
data_subset_gene$sample_cell = data_input$sample_cell
data_subset_gene$Treatment    =  data_input$Treatment                 


dir.create(paste0(filepath_cluster,'/Plots/Gene/'))
pdf(paste0(filepath_cluster,'/Plots/Gene/','violin_',gene,'_',celltype,'.pdf'), width =12)
plot = ggplot(data_subset_gene) + geom_violin(aes(x=label,y=gene,fill=label)) + theme_bw() +
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) + 
  theme(axis.text=element_text(angle=65,hjust=1))
print(plot)                                               
dev.off()

pdf(paste0(filepath_cluster,'/Plots/Gene/','boxplot_',gene,'_',celltype,'.pdf'), width =12)
plot = ggplot(data_subset_gene) + geom_boxplot(aes(x=label,y=gene,fill=label)) + theme_bw() +
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) + 
  theme(axis.text=element_text(angle=65,hjust=1))
print(plot)                                               
dev.off()

