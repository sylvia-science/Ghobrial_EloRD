library(wesanderson)
library(RColorBrewer)
library(scico)
library(pals)
library(nord)
library(palettetown)

folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/',
                'Plots/Paper/')


num_cluster = length(unique(Idents(data_harmony_run_label_remove)))
values = colorRampPalette(brewer.pal(12, "Accent"))(num_cluster)
values = colorRampPalette(wes_palette("Zissou1"))(num_cluster)

values = scico(num_cluster, palette = 'roma')
values = as.vector(ocean.delta(num_cluster))
values = rainbow(num_cluster) # Bright but not horrible
values = colorRampPalette(brewer.pal(11, "Paired"))(num_cluster) # Nice pastel feel
values = colorRampPalette(brewer.pal(8, "Dark2"))(num_cluster) # Too spooky


str = '_nolabel'
str = ''
labelTF = T
values  = hsv(seq(0, 1 - 1/num_cluster,length.out = num_cluster), .8, .85)
folder = paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/',
                'Plots/Paper/')

dir.create(folder, recursive = T)
pathName <- paste0(folder,'ClusterUmapAll_','hsv',str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_harmony_run_label_remove,pt.size = 0.7, reduction = "umap",label = labelTF, 
               cols =values)

plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()


values = colorRampPalette(brewer.pal(11, "Paired"))(num_cluster) # Nice pastel feel
pathName <- paste0(folder,'ClusterUmapAll_','Paired',str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_harmony_run_label_remove,pt.size = 0.7, reduction = "umap",label = labelTF, 
               cols =values)

plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()

###############
# T Cell
###############

str = '_nolabel'
str = ''
labelTF = T

for (i in 1:40){
  
  pokedex(10*i, 10)
}

Ident_order = c('Naive CD8+ T-cell', 'Naive CD4+ T-cell','Intermediate CD4+ T-cell','TSCM','Stim Naive CD4+ T-cell',
                'cTreg','eTreg','CD4+ TCM','TRM','Th2','Th17','aTh17','CCL5+ CD4+ T-cell',
                'CD8+ TCM','GZMK+ CD8+ T-cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell','GZMH+ GZMB+ CD8+ T-cell','TEMRA')
all(Ident_order %in% unique(Idents(data_run_subset_label) ))

Idents(data_run_subset_label) = factor(as.character(Idents(data_run_subset_label)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label)))


color = 'swablue'
values1 = colorRampPalette(ichooseyou(pokemon = color, spread = NULL))(15)
#values1 = values1[!(values1 %in% c('#F8F8F8','#E37E70'))]

color1 = 'swablu'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(5)

color2 = 'mewtwo'
values2 = colorRampPalette(ichooseyou(pokemon = color2, spread = NULL))(8)

values3 = colorRampPalette(ichooseyou(pokemon = 'magcargo', spread = NULL))(5)
values = c(values1,values2,values3)

pathName <- paste0(folder,'ClusterUmap_Tcell_',color1,'_',color2,str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_run_subset_label,pt.size = 0.5, reduction = "umap",label = labelTF, 
               cols =values)
plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()

################################
### Monocytes

str = '_nolabel'
str = ''
labelTF = T

unique(Idents(data_run_subset_label_remove))
       
Ident_order = c('cDC1','cDC2','prDC', 'sDC',
                'SELL+ CD14+ Mono','sMono','MIP1a+ CD14+ Mono','TGFb1+ CD14+ Mono',
                'sCD14+ Mono','IFN+ Mono','CD14+ CD16+ Mono','CD16+ Mono')
Ident_order[!(Ident_order %in% unique(Idents(data_run_subset_label_remove) ))]

#tmp = unique(Idents(data_run_subset_label_remove))[3]
Idents(data_run_subset_label_remove) = factor(as.character(Idents(data_run_subset_label_remove)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label_remove)))


color1 = 'venonat'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(4)

color2 = 'smoochum'
values2 = colorRampPalette(ichooseyou(pokemon = color2, spread = NULL))(8)
#values2 = values2[!(values2 %in% c('#F8F8F8'))]

values = c(values1,values2)

#color1 = 'magcargo'
#values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(num_cluster)
#values = values1

pathName <- paste0(folder,'ClusterUmap_Mono_',color1,'_',color2,str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_run_subset_label_remove,pt.size = 0.5, reduction = "umap",label = labelTF, 
               cols =values)
plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()

