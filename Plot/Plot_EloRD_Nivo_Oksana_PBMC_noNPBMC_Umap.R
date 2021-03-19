library(wesanderson)
library(RColorBrewer)
library(scico)
library(pals)
library(nord)
library(palettetown)

# Plot EloRD_Nivo_Oksana_PBMC_noNPBMC Umaps


## View palletes
for (i in 1:40){
  
  pokedex(10*i, 10)
}

########################
## Main clustering
########################

remove_list = c(0:100,'Remove')

data_harmony_run_label_remove = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
Idents(data_harmony_run_label_remove) = as.character(Idents(data_harmony_run_label_remove))

data_harmony_run_label_remove = data_harmony_run_label_remove[, !(Idents(data_harmony_run_label_remove) == 'T-cell' &
                                    data_harmony_run_label_remove@reductions[["umap"]]@cell.embeddings[,1] < 1)]
data_harmony_run_label_remove = data_harmony_run_label_remove[, !(Idents(data_harmony_run_label_remove) == 'CD14+ Mono' &
                                                                    data_harmony_run_label_remove@reductions[["umap"]]@cell.embeddings[,1] > -5)]

plot = DimPlot(data_harmony_run_label_remove,pt.size =0.3, reduction = "umap",label = TRUE,label.size = 3)
print(plot)

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
folder = paste0('/disk2/Projects/EloRD/Paper/Umap/')

dir.create(folder, recursive = T)
pathName <- paste0(folder,'ClusterUmapAll_',folder_name,'_','hsv',str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_harmony_run_label_remove,pt.size =0.3, reduction = "umap",label = labelTF, 
               cols =values, shape.by = NULL)

plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()


values = colorRampPalette(brewer.pal(11, "Paired"))(num_cluster) # Nice pastel feel

pathName <- paste0(folder,'ClusterUmapAll',folder_name,'_','Paired',str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 12
plot = DimPlot(data_harmony_run_label_remove,pt.size =0.3, reduction = "umap",label = labelTF, 
               cols =values,shape.by = NULL, raster = F)

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



Ident_order = c('Naïve CD8+ T-cell', 'Naïve CD4+ T-cell', 'aNaïve CD4+ T-cell','IFN+ T-cell','TSCM',
                'cTreg','eTreg','cCD4+ T-cell','CD4+ TCM','aCD4+ TCM','TRM','aTRM','Th2','Th17',
                'CD8+ TCM','GZMK+ CD8+ T-cell','GZMK+ GZMH+ CD8+ T-cell','CCL3+ CCL4+ GZMK+ CD8+ T-cell','GZMB+ GZMH+ CD8+ T-cell','TEMRA')
all(Ident_order %in% unique(Idents(data_run_subset_label) ))
data_idents = as.character(unique(Idents(data_run_subset_label)))
Ident_order[!(Ident_order %in%data_idents )]
data_idents[!(data_idents %in%Ident_order )]

data_run_subset_label_remove = data_run_subset_label[,as.character(Idents(data_run_subset_label)) %in% Ident_order]
plot = DimPlot(data_run_subset_label_remove,pt.size =0.3, reduction = "umap",label = TRUE,label.size = 3)
print(plot)

Idents(data_run_subset_label_remove) = factor(as.character(Idents(data_run_subset_label_remove)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label_remove)))

#data_run_subset_label_remove = data_run_subset_label_remove[, !(Idents(data_run_subset_label_remove) == 'T-cell' &
#                                                                  data_run_subset_label_remove@reductions[["umap"]]@cell.embeddings[,2] < 2.5)]
#data_run_subset_label_remove = data_run_subset_label_remove[, !(Idents(data_run_subset_label_remove) == 'CD8+ T-cell' &
#                                                                  data_run_subset_label_remove@reductions[["umap"]]@cell.embeddings[,2] < 2.5)]

color1 = 'kingdra'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(5)

color2 = 'mewtwo'
values2 = colorRampPalette(ichooseyou(pokemon = color2, spread = NULL))(9)

values2 [values2 == '#F8F8F8'] = '#892be0'
values3 = colorRampPalette(ichooseyou(pokemon = 'magcargo', spread = NULL))(6)
values = c(values1,values2,values3)

pathName <- paste0(folder,'ClusterUmap_Tcell_',folder_name,'_',color1,'_',color2,str,'.pdf')
pdf(file=pathName, width = 8,height = 8)

fontSize = 8
plot = DimPlot(data_run_subset_label_remove,pt.size = 0.3, reduction = "umap",label = labelTF, 
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
#################################
str = '_nolabel'
str = ''
labelTF = T

Ident_order = c('cDC1','cDC2','Neutrophil','MDPC',
                'HLA-DR-low SELL+ CD14+ Mono','SELL+ CD14+ Mono',
                'RNASE2-high SELL+ CD14+ Mono','TGFb1+ CD14+ Mono','CCL3+ CD14+ Mono',
                'IFN+ CD14+ Mono','CD14+ CD16+ Mono','CD16+ Mono','IFN+ CD16+ Mono')

all(Ident_order %in% unique(Idents(data_run_subset_label) ))
data_idents = unique(Idents(data_run_subset_label))
Ident_order[!(Ident_order %in%data_idents )]
data_idents[!(data_idents %in%Ident_order )]

Ident_order[!(Ident_order %in% unique(Idents(data_run_subset_label_remove) ))]

data_run_subset_label_remove = data_run_subset_label[,as.character(Idents(data_run_subset_label)) %in% Ident_order]
plot = DimPlot(data_run_subset_label_remove,pt.size =0.3, reduction = "umap",label = TRUE,label.size = 3)
print(plot)

data_run_subset_label_remove = data_run_subset_label_remove[,( data_run_subset_label_remove@reductions[["umap"]]@cell.embeddings[,2] >-4)]
data_run_subset_label_remove = data_run_subset_label_remove[,( data_run_subset_label_remove@reductions[["umap"]]@cell.embeddings[,2] <10)]



#tmp = unique(Idents(data_run_subset_label_remove))[3]
Idents(data_run_subset_label_remove) = factor(as.character(Idents(data_run_subset_label_remove)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label_remove)))


color1 = 'cascoon'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(4)


color2 = 'girafarig'
values2 = colorRampPalette(ichooseyou(pokemon = color2, spread = NULL))(9)

values = c(values1,values2)

#color1 = 'magcargo'
#values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(num_cluster)
#values = values1

pathName <- paste0(folder,'ClusterUmap_Mono_',folder_name,'_',color1,'_',color2,str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 8
plot = DimPlot(data_run_subset_label_remove,pt.size = 0.5, reduction = "umap",label = labelTF, 
               cols =values)
plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()

#########################
## NK
#########################
str = '_EloRD_PBMC'
str = ''
labelTF = T
data_idents =as.character( unique(Idents(data_run_subset_label)))

Ident_order = as.character(c('CD56br NK','aCD56br NK',data_idents[11],'CCL3+ NK',data_idents[1],data_idents[5],'IFN+ NK'))

all(Ident_order %in% as.character(unique(Idents(data_run_subset_label) )))
Ident_order[!(Ident_order %in%data_idents )]
data_idents[!(data_idents %in%Ident_order )]

Ident_order[!(Ident_order %in% unique(Idents(data_run_subset_label_remove) ))]

data_run_subset_label_remove = data_run_subset_label[,as.character(Idents(data_run_subset_label)) %in% Ident_order]
data_run_subset_label_remove = data_run_subset_label_remove[,( data_run_subset_label_remove@reductions[["umap"]]@cell.embeddings[,1] <5)]

plot = DimPlot(data_run_subset_label_remove,pt.size =0.3, reduction = "umap",label = TRUE,label.size = 3)
print(plot)


#tmp = unique(Idents(data_run_subset_label_remove))[3]
Idents(data_run_subset_label_remove) = factor(as.character(Idents(data_run_subset_label_remove)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label_remove)))

length(Ident_order)
color1 = 'rapidash'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(7)


values = c(values1)

#color1 = 'magcargo'
#values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(num_cluster)
#values = values1

pathName <- paste0(folder,'ClusterUmap_NK_',folder_name,'_',color1,str,'.pdf')
pdf(file=pathName, width = 8,height = 6)

fontSize = 8
plot = DimPlot(data_run_subset_label_remove,pt.size = 0.5, reduction = "umap",label = labelTF, 
               cols =values)
plot = plot + theme(
  legend.title = element_text( size = fontSize),
  legend.text = element_text( size = fontSize))
plot = plot +theme(axis.text=element_text(size=fontSize),
                   axis.title=element_text(size=fontSize,face="bold"))
print(plot)
dev.off()

##################
## HSC + B cell
##################

Ident_order = c('HSC','MEP','MP','GMP','DP',
                'Pro-B-cell','Pre-B-cell','Immature B-cell',
                'Naïve B-cell','Class-switched Memory B-cell','Plasma cell')



str = ''
labelTF = T
data_idents =as.character( unique(Idents(data_run_subset_label)))


Ident_order[!(Ident_order %in% data_idents)]
all(Ident_order %in% unique(Idents(data_run_subset_label) ))
Ident_order[!(Ident_order %in% data_idents )]
data_idents[!(data_idents %in% Ident_order )]

#Ident_order[!(Ident_order %in% unique(Idents(data_run_subset_label_remove) ))]

data_run_subset_label_remove = data_run_subset_label[,as.character(Idents(data_run_subset_label)) %in% Ident_order]
plot = DimPlot(data_run_subset_label_remove,pt.size =0.3, reduction = "umap",label = TRUE,label.size = 3)
print(plot)


#tmp = unique(Idents(data_run_subset_label_remove))[3]
Idents(data_run_subset_label_remove) = factor(as.character(Idents(data_run_subset_label_remove)),Ident_order)
num_cluster = length(unique(Idents(data_run_subset_label_remove)))


color1 = 'caterpie'
values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(5)
values1[ values1 == '#F8F8F8'] = '#48f0be'

color2 = 'starmie'
values2 = colorRampPalette(ichooseyou(pokemon = color2, spread = NULL))(6)

values = c(values1,values2)

#color1 = 'magcargo'
#values1 = colorRampPalette(ichooseyou(pokemon = color1, spread = NULL))(num_cluster)
#values = values1

pathName <- paste0(folder,'ClusterUmap_HSC_BCell_',folder_name,'_',color1,str,'.pdf')
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

