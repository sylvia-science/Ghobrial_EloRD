
data_input = data_run_subset_label
data_input = data_input[,!(Idents(data_input) %in% 
                             c('Erythrocyte',1:50,'Remove'))]

Idents(data_input) = as.character(Idents(data_input))


gene_list = read.csv('/home/sujwary/Desktop/scRNA/Data/GSEA/HALLMARK_INTERFERON_ALPHA_RESPONSE',header = F)
gene_list = gene_list$V1 
gene_list = as.character(gene_list[gene_list %in% rownames(data_run_subset)])

gene_list = list(c(gene_list))

data_run_subset_label = AddModuleScore(data_run_subset_label,gene_list, name = 'INTERFERON_ALPHA_RESPONSE')

val = cbind(as.character(Idents(data_run_subset_label)),as.numeric(data_run_subset_label$INTERFERON_ALPHA_RESPONSE1))
val = as.data.frame(val)
colnames(val) = c('label','INTERFERON_ALPHA_RESPONSE')
val$INTERFERON_ALPHA_RESPONSE = as.numeric(as.character(val$INTERFERON_ALPHA_RESPONSE))
aggregate(val$INTERFERON_ALPHA_RESPONSE, list(val$label), mean)

print(FeaturePlot(data_run_subset_label,pt.size = 0.7, features = c("INTERFERON_ALPHA_RESPONSE1")))
FontSize = 12

val = val[!(val$label %in% c(23, 25, 'Erythrocyte')),]

pathName =  paste0(filepath_cluster,'Plots/BoxPlot/','BoxPlot_Interferon','.png')
png(file=pathName,width=2400, height=600)

plot = ggplot(val, aes(x = label, y = INTERFERON_ALPHA_RESPONSE ,fill=label)) + 
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

###
list = c('Stim Naive CD4+ T-cell','Intermediate CD4+ T-cell','Naive CD4+ T-cell',
         'TSCM','CD4+ TCM','eTreg','Th2','Th17','aTh17')


gene_list = read.csv('/home/sujwary/Desktop/scRNA/Data/GSEA/HALLMARK_INTERFERON_ALPHA_RESPONSE',header = F)
gene_list = gene_list$V1 
gene_list = as.character(gene_list[gene_list %in% rownames(data_run_subset)])

pathName =  paste0(filepath_cluster,'Plots/HeatMap/','HeatMap_Interferon','.png')

png(file=pathName,width=2000, height=40*length(gene_list),res = 100)

plot = DoHeatmap(data_input, features = gene_list,slot = 'data')
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()


gene_list = 'CD4'
pathName =  paste0(filepath_cluster,'Plots/HeatMap/','HeatMap_CD4','.png')

png(file=pathName,width=2000, height=600,res = 100)

plot = DoHeatmap(data_input, features = gene_list,slot = 'data')
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()


