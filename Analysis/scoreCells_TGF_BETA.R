
data_input = data_harmony_run_label
#data_input = data_input[,!(Idents(data_input) %in% 
#                             c('Erythrocyte',1:50,'Remove'))]

#Idents(data_input) = as.character(Idents(data_input))


gene_list = read.csv('/home/sujwary/Desktop/scRNA/Data/GSEA/KEGG_TGF_BETA_SIGNALING_PATHWAY.txt',header = F)
gene_list = gene_list$V1 
gene_list = as.character(gene_list[gene_list %in% rownames(data_input)])

gene_list = list(c(gene_list))

data_input = AddModuleScore(data_input,gene_list, name = 'TGF_BETA')

val = cbind(as.character(Idents(data_input)),as.numeric(data_input$TGF_BETA1))
val = as.data.frame(val)
colnames(val) = c('label','TGF_BETA')
val$TGF_BETA = as.numeric(as.character(val$TGF_BETA))
val$cell = colnames(data_input)
val$sample = data_input$sample

write.csv(val,paste0(filepath_cluster,'/Data/TGF_BETA.csv'))

aggregate(val$TGF_BETA, list(val$label), mean)

print(FeaturePlot(data_input,pt.size = 0.7, features = c("TGF_BETA1")))
FontSize = 12

#val = val[!(val$label %in% c(23, 25, 'Erythrocyte')),]

pathName =  paste0(filepath_cluster,'Plots/BoxPlot/','BoxPlot_TGF_BETA1','.png')
png(file=pathName,width=2400, height=600)

plot = ggplot(val, aes(x = label, y = TGF_BETA ,fill=label)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, 0.6))+
  xlab("Cell Type") + ylab(paste0('Interferon Score'))+
  theme_classic()
# Box plot with dot plot
#plot = plot + geom_jitter(shape=16, position=position_jitter(0.2), color="black", size=2)

plot = plot + theme(
  #plot.title = element_text(hjust = 0.5,color="black", size=24, face="bold.italic"),
  axis.title.x = element_text(color="black", size=FontSize ),
  axis.title.y = element_text(color="black", size=FontSize),
  axis.text= element_text(color="black", size=FontSize),
)

plot = plot + theme(plot.title = element_text(hjust = 0.5))
print(plot)
dev.off()

data_matrix = data_input@assays[["RNA"]]@data
data_matrix_TGFb = data_matrix[rownames(data_matrix) %in% gene_list[[1]],]

path = paste0(filepath_cluster,'Data/matrix_',
              '','TGFb','.tsv')
write.table(data_matrix_TGFb, 
            file=path, 
            quote=FALSE, sep='\t')
