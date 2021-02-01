# merge the samples

# select highly variable genes

# run PCA and then 

# run BKNN instead of findneighbors. 
library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(harmony)
library(ggplot2)
library(SoupX)
library(sc)
library(scater)
library(dplyr)
library(scran)
library(reshape2)
library(stringr)
library(UpSetR)
library(edgeR)
library(SingleCellExperiment)
library(tibble)


source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')

source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
source('/home/sujwary/Desktop/scRNA/Code/Visualization/PlotCellPhoneDB.R')
source('~/Desktop/scRNA/Code/Integrate All/LoadHarmonySubset.R')

filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]
#metaData = metaData[metaData$`Sample Type` == 'PBMC',]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]

metaData = metaData[metaData$Sample !='BatchF',] # BatchF is crashing at quickCluster in scran norm

sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x

# Original ELoRD data with no PBMCs
sample_type = 'Harmony_AllSamples_Sample_Kit'
folder_name = 'AllSamples'
harmony_groupby = 'Sample_Kit'
base = '/home/sujwary/Desktop/scRNA/Output/Harmony/'

# EloRD data including our PBMCs and PBMCs for HCA
# sample_type = 'Harmony_AllSamples_PBMC_Sample_Kit'
# folder_name = 'AllSamples_PBMC'
# harmony_groupby = 'Sample_Kit'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

# # 
# EloRD data + our PBMCs + PBMCs for HCA
# sample_type = 'Harmony_AllSamples_PBMC_NPBMC_Sample_Kit'
# folder_name = 'AllSamples_PBMC'
# harmony_groupby = 'Sample_Kit'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

#EloRD data + our PBMCs + PBMCs for HCA + PBMCS/BMs from HCL
# sample_type = 'AllSamples_PBMC_NPBMC_HCL_Sample_kit_tech'
# folder_name = 'AllSamples_PBMC_NPBMC_HCL'
# harmony_groupby = 'Sample_kit_tech'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

#sample_type = 'Harmony_PBMC_Sample_Kit'


#folder_name = 'PBMC'





PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 


#patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'

#folder_name = 'AllSamplesDownsample'
#sample_list = Samples_runs$Samples[Samples_runs$Folder== folder]
#sample_list = unlist(strsplit(sample_list, ",")) 
#sample_list = trimws(sample_list, which = c("both"), whitespace = "[ \t\r\n]")
sample_list = metaData$Sample

#sample_list = sample_list[c(1, 3:8)]


folder = paste0(base,folder_name,
                '/Batch_',harmony_groupby,'/')
dir.create(folder,recursive = T)

sample_name = sample_list[1]
print(sample_name)
folder_input = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/', sample_name , '/')
data_1 = loadRData(paste0(folder_input,sample_name,'.Robj'))

run = F

if (run){
  data_list = vector(mode = "list",length = length(sample_list))
  data_list_norm = vector(mode = "list",length = length(sample_list))
  #for (i in 1:nrow(sampleParam)){
  for (i in 1:length(sample_list)){
    sample_name = sample_list[i]
    print(sample_name)
    folder_input = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/', sample_name , '/')
    data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
    data_i$sample = sample_name

    #path = paste0('/disk2/Projects/EloRD/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
    #cellIdents = read.csv(path,sep = ',',row.names = 1)
    #cellIdents$x = paste0(cellIdents$x, ' S', i)
    #data_i$CellType = cellIdents
    #data_i = data_i[,data_i$nFeature_RNA > 200]
    data_i = load_emptyDrops(data_i)
    data_i$is_cell[is.na(data_i$is_cell) ] = T # Not sure why there's NA values
    #print(unique(data_i$is_cell))
    print(sum(!data_i$is_cell == T))
    #next
    # remove empty drops
    data_i = data_i[,data_i$is_cell]
    
    data_i = data_i[rownames(data_1),]
    # Remove MT > 15 already done
    if (!is.na(downsample)){
      downsample = sub("_.*", "", downsample)
      cellnames = colnames(data_i)
      cellnames = sub("_.*", "", cellnames)
      data_i = data_i[,cellnames %in% downsample]
      #browser()
    }
    
    print(ncol(data_i))
    if (ncol(data_i) > 100){
      data_list[[i]] = data_i
      data_list_norm[[i]] = ScranNorm(data_i)
      
    }
    
  }
  
  data_list_norm =data_list_norm[lengths(data_list_norm) != 0]
  data_merge = merge(x =  data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
  data_merge = addMetaData(data_merge, metaData)
  #data_merge = load_emptyDrops(data_merge)
  data_merge = load_Doublets(data_merge)
  #data_merge = load_CellLabel(data_merge)
  #data_merge$GeneralCellType = str_match(data_merge$CellType, "(^.+)\\s")[, 2]
  data_merge$kit = data_merge$`10X kit`
  data_merge$split_var = ''
  cell_names_all = sub("_.*", "", colnames(data_merge))
  data_merge$cell_sample = paste0(data_merge$sample ,' ',cell_names_all)
  data_merge$sample_cell = paste0(cell_names_all,' ',data_merge$sample )
   
  data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  #data_harmony_run = ScaleData(data_harmony_run, vars.to.regress = c('kit'))
  data_merge_run = ScaleData(data_merge_run)
  
  data_merge_run = RunPCA(data_merge_run,npcs = 70)
  
  reduction = 'pca'


  visualize_PCA(data_merge_run,folder,70,reduction)
  
  
  PCA_dim = 60
  resolution_val = 1
  data_merge_run = RunUMAP(reduction = "pca",data_merge_run, dims = 1:PCA_dim)
  data_merge_run = FindNeighbors(data_merge_run, reduction = "pca", dims = 1:PCA_dim)
  data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  
  data_merge_run$FeatureLessThan400 = data_merge_run$nFeature_RNA < 400
  data_merge_run$kit = data_merge_run$`10X kit`
  
  folder = paste0(base,folder_name,
                  '/Merge','/')
  path = paste0(folder,'data_merge_run','.Robj')
  save(data_merge_run,file= path)
  data_merge_run = loadRData((path))
  
  data_merge_run$split_var = ''
  
  metaData = read_excel(filename_metaData)
  data_merge_run = addMetaData(data_merge_run, metaData)
  data_merge_run$kit = data_merge_run$`10X kit`
  #
  groupBy_list = c('sample','Diagnosis','kit',
                   'Treatment','Batch','LowCount',
                   'FeatureLessThan400','Sample Type','Technology')
  featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
  splitBy_list = NA
  
  data = as.data.frame(data_merge_run@assays[["RNA"]]@counts)
  DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
  data_merge_run$DR = DR
  

  plotAll(data_merge_run, folder = folder,
          sample_name,sampleParam = NA,
          cell_features = cell_features, plot_PCA = F,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = F,
          groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
          PCA_dim = PCA_dim,resolution_val = resolution_val)
  
  
  harmony_dim = PCA_dim
  data_harmony_run = RunHarmony(data_merge_run,group.by.vars =  c("kit",'sample','Technology'),
                                dims.use = 1:harmony_dim)
  
  folder = paste0(base,folder_name,
                  '/Batch_',harmony_groupby,'/')
  reduction = 'harmony'
  visualize_PCA(data_harmony_run,folder,harmony_dim,reduction)
  
  data_harmony_run = RunUMAP(reduction = "harmony",data_harmony_run, dims = 1:harmony_dim)
  data_harmony_run = FindNeighbors(data_harmony_run, reduction = "harmony", dims = 1:harmony_dim)
  data_harmony_run = FindClusters(data_harmony_run,resolution = resolution_val)
  
  
  plot = DimPlot(data_harmony_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
  path = paste0(folder,'data_run','.Robj')
  save(data_harmony_run,file= path)
  

}else{
  resolution_val = 3
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  
  path = paste0(folder,'data_run','.Robj')
  data_harmony_run = loadRData(path)
  tmp = data_harmony_run@meta.data[paste0('RNA_snn_res.', resolution_val)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_harmony_run) = tmp
  data_harmony_run$orig.ident = Idents(data_harmony_run)
  
  cell_names_all = sub("_.*", "", colnames(data_harmony_run))
  data_harmony_run$cell_sample = paste0(cell_names_all,' ',data_harmony_run$sample )
  data_harmony_run$sample_cell = paste0(data_harmony_run$sample,' ',cell_names_all )
  
  cell_names_main = as.character(data_harmony_run$cell_sample)

  path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3//Data/Labels.csv'
  data_harmony_run = addOldLabels(path,data_harmony_run, 'OldCellType')
  #data_harmony_run_label = LoadSubset(data_harmony_run_label,sampleParam_combine, folder)
  data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)
  
  celltype = 'Mono_DC'
  resolution_val_subset = 1.6
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run_PC30','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  #data_run_subset_label = data_run_subset
  #Idents(data_run_subset_label) = paste0(celltype,' ', Idents(data_run_subset_label))
  
  
  Ident_main = colnames(data_harmony_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_harmony_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_harmony_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_harmony_run_label) = newIdents2
  ##
  
  celltype = 'NK_RemoveRiboFeatures'
  #celltype = 'NK'
  
  resolution_val_subset = 3
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run_PC30','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  Ident_main = colnames(data_harmony_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_harmony_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_harmony_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_harmony_run_label) = newIdents2
  
  
  #
  
  celltype = 'T Cell'
  resolution_val_subset = 3.5 # main run
  PCA_subset = 30
  #resolution_val_subset = 4 # main run
  #PCA_subset = 40
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run_PC',PCA_subset,'.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  #data_run_subset_label = data_run_subset
  #Idents(data_run_subset_label) = paste0(celltype,' ', Idents(data_run_subset_label))
  
  
  Ident_main = colnames(data_harmony_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_harmony_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_harmony_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_harmony_run_label) = newIdents2
  
  
  
  #(0,48)
  #data_harmony_run = FindClusters(data_harmony_run,resolution = resolution_val)
  
  data_harmony_run$sample_cell = paste0(data_harmony_run$Sample,'_', colnames(data_harmony_run))
  
  TSCM_res4 = read.csv(paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/T Cell/Cluster/PCA30/res3.5/'
                              ,'Data/','TSCM_res4','.csv'))
  
  new_Idents = as.character(Idents(data_harmony_run_label))
  new_Idents[colnames(data_harmony_run_label) %in% TSCM_res4$x] = 'TSCM'
  Idents(data_harmony_run_label) = new_Idents
  
}


data_harmony_run$GeneralCellType = ''
metaData = read_excel(filename_metaData)


table(data_harmony_run$sample)
table(data_harmony_run$Treatment)
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
dir.create(filepath_cluster,recursive = T)


#data_harmony_run = load_emptyDrops(data_harmony_run)
data_harmony_run = load_Doublets(data_harmony_run)
print(unique(data_harmony_run$GeneralCellType))
#data_harmony_run$GeneralCellType = str_match(data_harmony_run$CellType, "(^.+)\\s")[, 2]
#data_harmony_run$split_var = ''
data_harmony_run$FeatureLessThan400 = data_harmony_run$nFeature_RNA < 400
data_harmony_run = addMetaData(data_harmony_run, metaData)
data_harmony_run$kit = data_harmony_run$`10X kit`

##
data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)
##

remove_list = c('Remove','14','32','Erythrocyte','DC/T-Cell DBL',
                'Mono/CD8+ T Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL','42','Plasma Cell','41','',0:40,
                'Remove 29','Remove 30','CD14+ Mono/T-cell DBL','CD14+ Mono/CD8+ T-cell DBL','CD16+ Mono/T-cell DBL')

remove_list = c('Remove','14','32','Erythrocyte',
                '42','Plasma Cell','41','',0:40,
                'Remove 29','Remove 30')

remove_list = c('13','30','32','41','49','42','Remove')

data_harmony_run_label_remove = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
Idents(data_harmony_run_label_remove) = as.character(Idents(data_harmony_run_label_remove))
unique(Idents(data_harmony_run_label_remove))
#data_harmony_run_label_remove = data_harmony_run_label[,data_harmony_run_label$`Sample Type` == 'PBMC']


#newIdents = as.character(Idents(data_harmony_run))
#newIdents[newIdents == 15] = 'GZMB+GZMH+'
#newIdents[newIdents %in%  c(6,14)] = 'GZMK+'
#Idents(data_harmony_run_label) = newIdents

#data_harmony_run_label = data_harmony_run_label[,data_harmony_run_label$is_cell]
#data_harmony_run_label = data_harmony_run_label[,Idents(data_harmony_run_label) != 'Remove']


plot = DimPlot(data_harmony_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

plot = DimPlot(data_harmony_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

plot = DimPlot(data_harmony_run_label_remove,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_harmony_run_subset = data_harmony_run[,!(data_harmony_run$OldCellType == '')]
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy','OldCellType_label','.png'))
png(file=pathName,width=1000, height=1000)

plot = DimPlot(data_harmony_run_subset, pt.size = 0.7, reduction = "umap",label = T,
               group.by  = 'OldCellType', label.size = 8)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 24)
)
print(plot)
dev.off()


cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

groupBy_list = c('sample','Diagnosis','kit',
                 'Treatment','Batch','LowCount',
                 'Doublet','GeneralCellType','FeatureLessThan400','Sample Type')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score','DR')
splitBy_list = NA

#data = as.data.frame(data_harmony_run@assays[["RNA"]]@counts)
#DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
#data_harmony_run$DR = DR

plotAll(data_harmony_run, folder = folder, 
        sample_name,sampleParam = NA,
        cell_features = cell_features, plot_PCA = F,
        label_TF = F,integrate_TF = F,  DE_perm_TF = T, 
        clusterTF =F, markersTF = F,
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val)


plotAll(data_harmony_run_label, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,plot_PCA = F,
        label_TF = F,integrate_TF = T,  DE_perm_TF = F, 
        clusterTF =F, markersTF = T, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val,str = '_label')

data_harmony_run_label_remove

plotAll(data_harmony_run_label_remove, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,plot_PCA = F,
        label_TF = F,integrate_TF = T,  DE_perm_TF = F, 
        clusterTF =F, markersTF = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val,str = '_label_remove')

## Print old cell type for cluster
cell_table = table( data_harmony_run$OldCellType,Idents(data_harmony_run))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','OldCellTypeByIdent_Table.csv'))

cell_table = table( data_harmony_run$sample,Idents(data_harmony_run))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','sampleByIdent_Table.csv'))

cell_table = table( as.character(data_harmony_run$sample))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','NumCellPerSample.csv'))


for (cluster in sort(unique(Idents(data_harmony_run)))){
  print(cluster)
  data_harmony_run_subset = data_harmony_run[,Idents(data_harmony_run) == cluster]
  print(table(data_harmony_run_subset$OldCellType))
}

########################
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
#cell_features = cell_features[cell_features$Cell == 'monocyte_FCGR3A',]
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )

PlotKnownMarkers(data_harmony_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F)

cell_features = cell_features[cell_features$Plot_marker == 1,]
PlotKnownMarkers(data_harmony_run, folder = paste0(filepath_cluster,'Cell Type/AllGenes/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F, plotAll = T)
print(DimPlot(data_harmony_run, reduction = "umap", label = TRUE, pt.size = .1))

gene_list = c('IRF4','IRS1','CARM1','CFLAR','PRDM1','MARCH5','FURIN','HERPUD1')

FeaturePlot_GeneList(data_harmony_run,gene_list,
                     folder = paste0(filepath_cluster,'Cell Type/Mahshid/'),
                     FeaturePlotFix = T,str = '')

StackedVlnPlotHelper(data_harmony_run_label,gene_list,
                                folder_heatMap = paste0(filepath_cluster,'Plots/Violin/'),
                                filename = paste0('AXL_GAS6_AllCluster.png'))
  

tmp  = data_harmony_run_label$sample[Idents(data_harmony_run_label) == 'B Cell' & data_harmony_run_label$Treatment == 'baseline']
table(as.character(tmp ))

sum(Idents(data_harmony_run_label) == 'Plasma Cell' & data_harmony_run_label$Treatment == 'baseline')

group = 'CellType'
pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',30,'_res',resolution_val,'_GroupBy',group,'.png'))
png(file=pathName,width=5000, height=1500)
plot = DimPlot(data_harmony_run,pt.size = 1, reduction = "umap",label = FALSE,group.by  = group)

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20)
)


#plot = plot + scale_color_manual(values=color_list)

print(plot)

dev.off()



####################
### Plot DE
####################
str = ''

filepath = paste0(filepath_cluster
                  ,'Features',str,
                  '.csv')
DE = read.csv(file = filepath)
data_input = data_harmony_run


DE = DE[DE$p_val_adj < 0.05,]
DE = DE[!grepl("MT-", DE$gene),]
DE = DE[!grepl("^RP[SL]", DE$gene),]

cluster_list = unique(DE$cluster)

#cluster_list = cluster_list[4:length(cluster_list)]
for (cluster in cluster_list){
  DE_cluster = DE[DE$cluster == cluster,]
  ranks <- DE_cluster$avg_logFC
  names(ranks) <- DE_cluster$gene
  
  DE_cluster = DE_cluster[DE_cluster$avg_logFC > 0 ,]
  DE_cluster = DE_cluster[order(-DE_cluster$avg_logFC),]
  DE_cluster = DE_cluster[1:20,]
  
  DE_cluster = DE_cluster[rowSums(is.na(DE_cluster)) != ncol(DE_cluster),]
  gene_list = DE_cluster
  
  
  gene_list$Cell = cluster
  gene_list$Plot_marker = 1
  gene_list$Markers = gene_list$gene
  path = paste0(filepath_cluster,'DE/Cluster_',cluster,'VsAll',str,'/')
  PlotKnownMarkers(data_input, 
                   folder = path, 
                   cell_features = gene_list,
                   plotType ='FeaturePlotFix' , 
                   prefix_logFC = T, str = '',markerSize = 1)
  
  
  
}

#######################
## DE
#######################
# Detection rate
#data = as.data.frame(data_harmony_run_label@assays[["RNA"]]@counts)
#DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
#data_harmony_run_label$DR = DR
source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')

data_harmony_run_label = ScaleData(data_harmony_run_label, features = rownames(data_harmony_run_label))


celltype_list = unique(as.character(Idents(data_harmony_run_label)))
celltype_list

celltype_num = sort(table(as.character(Idents(data_harmony_run_label))))

celltype_list = c('cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')

celltype_list = c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell')


celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')



celltype_list = c('Cytotoxic NK')

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')

celltype_list = c('dMono','Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype_list = c('GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype = celltype_list[1]
celltype_list = c('Inhibitory NK')


#celltype_list = c('Plasma Cell')
#celltype_list=  c('CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')
for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0('baseline Yes ',celltype)
  ident2 = paste0('baseline No ',celltype)
  
  #ident1 = paste0('CD16+ Mono')
  #ident2 = paste0('CD14+ Mono')
  
  DE_input = data_harmony_run_label
  DE_input = renameCells(DE_input,idents = c('cDC1','cDC2'),newident = 'DC')
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  #DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
  #                          DE_input$Best_Overall_Response, ' ', Idents(DE_input))
  DE_input$DE_ident = paste0(DE_input$Treatment,' ',DE_input$'Dexa or not', ' ', Idents(DE_input))
  #DE_input$DE_ident = paste0(Idents(DE_input))
  #DE_input = DE_input[,DE_input$Treatment =='baseline']
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
  ncol(DE_input)
  unique(DE_input$DE_ident)
  
  if (ncol(DE_input) < 100){
    next
  }
  
  #DE_input_sce = as.SingleCellExperiment(DE_input)
  
  data = as.data.frame(DE_input@assays[["RNA"]]@counts)
  # detection rate:fraction of genes expressed in a cell
  DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
  DE_input$DR = DR 
  
  kit = factor(DE_input$kit)
  ident = factor(as.character(Idents(DE_input)))
  DE_input$ident = ident
  DE_ident = factor(DE_input$DE_ident)
  Patient = factor(DE_input$`Patient Number`)
  DE_input$Patient = Patient
  Treatment = factor(DE_input$Treatment)
  
  if (length(unique(DE_ident)) == 1){
    next
  }
  library(robustbase)
  library(DESeq2)
  
  if (length(unique(kit)) == 1){
    formula = ~  DE_ident + DR + Patient
  }else{
    formula = ~  DE_ident + DR + kit + Patient
  }
  #formula = ~  DE_ident + DR + kit
  design <- model.matrix(formula)
  colnames(design) <- gsub("DE_ident", "", colnames(design))
  colnames(design)
  is.fullrank(design)
  design_fr = fullRank(design)
  is.fullrank(design_fr)
  colnames_old = colnames(design)
  colnames_fr = colnames(design_fr)
  
  colnames_old[!(colnames_old %in% colnames_fr)]
  
  keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                       min.count = 1,min.total.count=10, 
                       large.n = 10,min.prop = 0.1)
  subfolder = paste0(ident1,' Vs ',ident2)
  
  result_DESeq2 = runDESeq2(DE_input,design_fr,contrast, keep, 
                     folder_output = filepath_cluster, subfolder = subfolder)
  
  dds = result_DESeq2[[1]]
  res_DESeq2 = result_DESeq2[[2]]
  #plotDispEsts(dds)
  
  res_DESeq2 = res_DESeq2[order(res_DESeq2$log2FoldChange),]
  res_DESeq2 = res_DESeq2[res_DESeq2$pvalue < 0.05,]
  
  res_DESeq2 = res_DESeq2[!grepl("MT-", res_DESeq2$gene),]
  res_DESeq2 = res_DESeq2[!grepl(" ?RP\\w+ ?", res_DESeq2$gene),]
  
  folder = paste0(filepath_cluster,'/DE/DESeq2/',subfolder,'/')
  dir.create(folder,recursive = T)
  pathName <- paste0(folder,paste0('HeatMap','.png'))
  #png(file=pathName, height=12, width=20)
  #png(file=pathName)
  len = nrow(res_DESeq2)
  if (nrow(res_DESeq2) > 200){
    res_DESeq2 = res_DESeq2[c(1:100,(len - 100):len),]
    res_DESeq2 = res_DESeq2[rowSums(is.na(res_DESeq2)) != ncol(res_DESeq2),]
    res_DESeq2 = unique(res_DESeq2)
  }
  
  
  png(file=pathName,width=2000, height=40*nrow(res_DESeq2),res = 100)
  
  plot = DoHeatmap(DE_input, features = res_DESeq2$gene, group.by = 'DE_ident')
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20))
  print(plot)
  dev.off()

}
########################



#######################
sum(res_DESeq2$padj < 0.05, na.rm = T)
sum(res_DESeq2$padj > 0.05, na.rm = T)
nrow(res_DESeq2)

res_DESeq2_sig = res_DESeq2$gene[res_DESeq2$padj < 0.05 & !is.na(res_DESeq2$padj)]


  
res <- results(dds, 
               name = names[2],
               alpha = 0.05,  pAdjustMethod	 ='BH')
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_sig = res_tbl[res_tbl$padj < 0.05 & !(is.na(res_tbl$padj)),]
res_tbl_sig


Idents(DE_input) = DE_input$DE_ident
Features = FindMarkers(DE_input, ident.1 = ident1, ident.2 = ident2
                       ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, slot = "counts")

Features_sig = Features[Features$p_val_adj< 0.05,]

intersect(res_tbl_sig$gene, rownames(Features_sig))
#  Identifies differentially expressed genes between two groups of cells 
#  using a hurdle model tailored to scRNA-seq data. 
# Utilizes the MAST package to run the DE testing.

# MAST
markers = FindMarkers(data_harmony_run_label,  
                      ident.1 = "CD8+ T Cell", ident.2 = "T Cell", 
                      latent.vars = 'kit', test.use = "MAST")

##
# edgeR
data_harmony_run_label_sce = as.SingleCellExperiment(data_harmony_run_label)
## Convert to DGEList, calculate logCPMs
dge <- scran::convertTo(data_harmony_run_label_sce, type = "edgeR")
plotMDS(dge)

design <- model.matrix(~kit+DR)

dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)
topTags(lrt)

##
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)



# voom limma
dge <- DGEList(data_harmony_run_label@assays[["RNA"]]@counts, 
               group = data_harmony_run_label$kit)
dge <- calcNormFactors(dge)

design <- model.matrix(~ kit + DR)

vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

hist(tt$P.Value, 50)
hist(tt$adj.P.Val, 50)
#limma::plotMDS(dge, col = as.numeric(as.factor(data_harmony_run_label$kit)), pch = 19)
plotMD(fit)

# Cluster comparisions using design matrix made from cell idetents
design <- model.matrix(~ 0 + ident, data = colData(data_harmony_run_label_sce))
colnames(design) <- gsub("ident", "", colnames(design))
colnames(design)

y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)

## Perform pairwise comparisons
nclust <- length(unique(data_harmony_run_label_sce$ident))
all.results <- all.pairs <- list()
counter <- 1

for (i in seq_len(nclust)) {
  for (j in seq_len(i - 1L)) {
    con <- integer(ncol(design))
    con[i] <- 1
    con[j] <- -1
    fit2 <- limma::contrasts.fit(fit, con)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    
    res <- limma::topTable(fit2, number = Inf, sort.by = "none")
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(i, j)]
    counter <- counter + 1L
    
    ## Also filling the reverse comparison.
    res$logFC <- -res$logFC
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(j, i)]
    counter <- counter + 1L
  }
}

## Combine results across all pairwise tests
all.pairs <- do.call(rbind, all.pairs)
combined <- scran::combineMarkers(all.results, all.pairs, 
                                  pval.field = "P.Value",
                                  pval.type = "any")
head(combined[["cluster1"]])

markers <- scran:::findMarkers(data_harmony_run_label_sce, design=design)

#scran:::findMarkers()

#####################
## Heatmap
######################


data_harmony_run_label = ScaleData(data_harmony_run_label,features =  rownames(data_harmony_run_label))

remove_list= c('HSC','Pro Erythrocyte','MDPC','MK','Erythrocyte','Pro B Cell','Pre B Cell','Remove',
  'CMPC','42','GMPC','Plasma Cell','B Cell')

data_harmony_run_label = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]

gene_list = gene_list = c('AXL','GAS6')
folder = paste0(filepath_cluster,'/Plots/Heatmap/')
dir.create(folder,recursive = T)
pathName <- paste0(folder,paste0('HeatMap', '_AXL_GAS6','.pdf'))
pdf(file=pathName, height=12, width=20)
#png(file=pathName)
#png(file=pathName,width=2000, height=1000,res = 100)

plot = DoHeatmap(data_harmony_run_label, features = gene_list)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()

######################
cluster_list = levels(unique(Idents(data_harmony_run)))

for (cluster in cluster_list){
  data_subset = data_harmony_run[,Idents(data_harmony_run) == cluster]
  print(cluster)
  cell_percent =  100*table(data_subset$GeneralCellType)/ncol(data_subset)
  print(cell_percent)
  print(max(cell_percent))
  print('')
}

for (cluster in cluster_list){
  data_subset = data_harmony_run[,Idents(data_harmony_run) == cluster]
  print(cluster)
  sample_percent =  100*table(data_subset$sample)/ncol(data_subset)
  print(sample_percent)
  print(max(sample_percent))
  print('')
}

filepath_cluster = paste0( folder, 'Cluster/', 'PCA',30,'/res',resolution_val,'/' )
k_num = 30
dim_num = 30
folder_output = paste0(filepath_cluster,'Entropy/')
dir.create(folder_output)
folder_output = paste0(folder_output,'Harmony')


#data_harmony_run_label_downsample = data_harmony_run_label[,cellNames]

compute_entropy_Seurat(data_harmony_run, 
                       corrected_assay = 'harmony',folder = folder_output,
                       k_num = k_num, dim_num = dim_num)


file_entropy = paste0(folder_output,'_k',k_num ,'_entropy.csv')
entropy = read.csv(file_entropy,sep = ',')

entropy$Method = 'Harmony'

entropy$batch_entropy <- NULL
plotEntropy(entropy,folder_output)


##########################
sample_summary =table(data_harmony_run$sample)
write.csv(sample_summary,file = paste0(filepath_cluster,'Stats/','celltype_summary','.csv'))

cellPerSample =table(data_harmony_run_label$sample, Idents(data_harmony_run_label))
write.csv(cellPerSample,file = paste0(filepath_cluster,'Stats/','cellPerSample','.csv'))
####################
## Doublets
####################
doubletMethod = 'Doublet3Methods'
mean_doublets = 100*sum(data_harmony_run$Doublet3Methods)/ncol(data_harmony_run)
title_str = paste0('% Doublets: ',specify_decimal(mean_doublets,2))
pathName <- paste0(filepath_cluster,paste0(doubletMethod,'_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_harmony_run,
              pt.size = ifelse(data_harmony_run$Doublet3Methods == T, 1, 0.5), cols = c('navajowhite2','tomato2'), 
              reduction = "umap",label = F, group.by  = doubletMethod) + 
              ggtitle(title_str)
plot = plot +  DimPlot(data_harmony_run,pt.size = 0,
               reduction = "umap",label = T, group.by  = 'GeneralCellType')
print(plot)
dev.off()


doublet_data = data.frame(cbind(data_harmony_run$scrublet,
                                data_harmony_run$doublet_finder, 
                                data_harmony_run$scran_doublet, 
                                data_harmony_run$scds_doublet))
colnames(doublet_data) = c('scrublet','doublet_finder','scran_doublet','scds_doublet')

doublet_data = doublet_data*1

pathName <- paste0(filepath_cluster,paste0('Stats/','doublet_upset','.pdf'))
pdf(file=pathName)
plot = upset(doublet_data, sets =  c('scrublet','doublet_finder','scran_doublet','scds_doublet'),
      point.size = 3.5, line.size = 2, order.by = "freq",
      mainbar.y.label = "Doublet Intersections", sets.x.label = "Doublets Per method")
print(plot)
dev.off()

##########################
##
##########################
pathName <- paste0(filepath_cluster,paste0('Empty_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_harmony_run,
               pt.size = ifelse(data_harmony_run$is_cell == T, 0.5, 1), cols = c('tomato2','navajowhite2'), 
               reduction = "umap",label = F, group.by  = 'is_cell')
plot = plot +  DimPlot(data_harmony_run,pt.size = 0,
                       reduction = "umap",label = T, group.by  = 'GeneralCellType')
print(plot)
dev.off()

downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
# Plot scatterplots of gene


gene_list = c('HBB','IGKC','IGLL5')

gene1 = 'HBB'
for (gene1 in gene_list){
  plotGeneScatter(data_harmony_run,gene1, 'PTPRC')
  plotGeneScatter(data_harmony_run,gene1, 'CD3D')
  plotGeneScatter(data_harmony_run,gene1, 'CD14')
  plotGeneScatter(data_harmony_run,gene1, 'FCGR3A')
  plotGeneScatter(data_harmony_run,gene1, 'NKG7')
}


## Get T cell, mono stats
path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/Mono_DC/Cluster/PCA30/res1.6/Stats/TcellMono_cells'
TcellMono_cells = read.csv(path)

unique(TcellMono_cells$samples)

orig_cell_labels = as.character(Idents(data_harmony_run_label))

orig_cell_labels[data_harmony_run_label$cell_sample %in% TcellMono_cells$cell_sample] = 'T_Mono'
orig_cell_labels[orig_cell_labels == 'CD8+ T Cell'] = 'T Cell'
orig_cell_labels[orig_cell_labels == 'CD14+ Mono' | orig_cell_labels == 'CD16+ Mono'] = 'Mono'

df = data.frame(cbind(orig_cell_labels,data_harmony_run_label$sample))

colnames(df) = c('Cell','Sample')

df = df[df$Cell %in% c('T Cell','Mono','T_Mono'),]
df$Cell = as.character(df$Cell)
df$Sample = as.character(df$Sample)
#df  = df[df$Sample %in% unique(TcellMono_cells$samples),]
df_stats = t(table(df))

write.csv(df_stats,file = paste0(filepath_cluster,'Stats/','T_Mono_Stats.csv'))


#######################################3
cell_list = c("Naive CD8+ T-cell","cTreg" ,"GZMH+ GZMB+ CD8+ T-cell","CD8+ TCM",
              "GZMK+ CCL3+ CCL4+ CD8+ T-cell", "Intermediate CD4+ T-cell" ,
              "TSCM"    ,"Naive CD4+ T-cell"  ,"eTreg"  , "CD4+ TCM" ,"Th2",  "TEMRA" ,                   
              "GZMK+ CD8+ T-cell" ,"Stim Naive CD4+ T-cell" , "CCL5+ CD4+ T-cell"  ,"aTh17" , "Th17","TRM",
              'pDC', "cDC2" , "sDC"   ,                  
              "prDC"  ,  "cDC1"  , 
              "CD16+ Mono" ,               "TGFb1+ CD14+ Mono",      
              "SELL+ CD14+ Mono" ,          "CD14+ CD16+ Mono"   ,     
              "IFN+ Mono"      ,           "MIP1a+ CD14+ Mono",
              "Tgd"  ,"aCX3CR1+ GZMB+ CD56dim" ,   "aCXCR4+ CD56dim"    
              ,"aGZMK+ CCL3+ CD56dim",  "cCD56dim" ,"CD56bright"   ,"NFkB-high")  

unique(Idents(data_harmony_run_label))

cell_list[!(cell_list %in% unique(Idents(data_harmony_run_label)))]


data_input = data_harmony_run_label[,Idents(data_harmony_run_label) %in% cell_list]

data_input = FindVariableFeatures(data_input, selection.method = "vst", nfeatures = 2000)

path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit_',
              'TCell_NK_Mono','.tsv')
data_matrix = data_input@assays[["RNA"]]@data
data_matrix = data_matrix[rownames(data_matrix) %in% data_input@assays[["RNA"]]@var.features,]
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')

path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit_',
              'TCell_NK_Mono_cellTypes','.tsv')
celltypes = as.data.frame(as.character(Idents(data_input)))
colnames(celltypes) = 'cellType'
celltypes$sample_cell = data_input$sample_cell
write.table(celltypes, 
            file=path, 
            quote=FALSE, sep='\t')

sort(table(celltypes$cellType))

## Rerun data and Save
harmony_dim = 30
data_input = RunUMAP(reduction = "harmony",data_input, dims = 1:harmony_dim)
data_input = FindNeighbors(data_input, reduction = "harmony", dims = 1:harmony_dim)
data_input = FindClusters(data_input,resolution = 3)

data_input$ident_orig = as.character(Idents(data_harmony_run_label)[Idents(data_harmony_run_label) %in% cell_list])
plot = DimPlot(data_input,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

plot = DimPlot(data_input,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8,group.by  = 'ident_orig')
print(plot)

path = paste0(folder,'data_T_NK_Mono_run','.Robj')
save(data_input,file= path)

sample = as.character(data_input$sample)
ident = as.character(Idents(data_input))
UmapCoord = data_input@reductions[["umap"]]@cell.embeddings
ident_orig = Idents(data_harmony_run_label)[Idents(data_harmony_run_label) %in% cell_list]

output = as.data.frame(sample)
output$ident = ident
output$ident_orig = data_input$ident_orig

output = cbind(output,UmapCoord)
#output = data_harmony_run_label@meta.data

dir.create(paste0(filepath_cluster,'data/'))
path = paste0(filepath_cluster,'data/TCell_NK_Mono_cellTypes','.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)


##########33333333333333333333333333
path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit',
              '.tsv')
data_matrix = data_harmony_run@assays[["RNA"]]@data
data_matrix = data_matrix[rownames(data_matrix) %in% data_harmony_run@assays[["RNA"]]@var.features,]
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')

# Save all data
#data_harmony_run_label = data_harmony_run_label[,Idents(data_harmony_run_label) %in% c('NK','CD14+ Mono')]
data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)

data_mere_run_label = data_main_label

data_harmony_run_label_input = data_harmony_run_label[,data_harmony_run_label$Treatment == 'baseline']
library("org.Hs.eg.db") # remember to install it if you don't have it already

data_matrix = as.data.frame(data_harmony_run_label_input@assays[["RNA"]]@data)
gene_ensemble= mapIds(org.Hs.eg.db, keys = rownames(data_matrix), keytype = "SYMBOL", column="ENSEMBL", multiVals = "asNA")
gene_ensemble = gene_ensemble[!is.na(gene_ensemble)]
gene_ensemble = gene_ensemble[!duplicated(unname(gene_ensemble)) ] # Keeps only first dup
data_matrix = data_matrix[names(gene_ensemble),]
rownames(data_matrix) = unname(gene_ensemble)

path = paste0(filepath_cluster,'/data/')
dir.create(path,recursive = T)

str = '_baseline'

write.csv(data_matrix, 
            file=paste0(path,'Harmony_AllSamples_Sample_Kit',str,'.csv'), 
            quote=FALSE, row.names = T)

write.table(data_matrix, 
            file=paste0(path,'Harmony_AllSamples_Sample_Kit',str,'.tsv'), 
            quote=FALSE, sep='\t', row.names = T, col.names=NA)

path = paste0(filepath_cluster,'/data/')
labels = as.character(Idents(data_harmony_run_label))
write.csv(labels, file = paste0(path,'Labels.csv'), )

labels = as.data.frame(colnames(data_harmony_run_label))
colnames(labels) = 'Cell'
labels$sample = (as.character(data_harmony_run_label$sample))

labels$cell_type = (as.character(Idents(data_harmony_run_label)))
labels$doublet = data_harmony_run_label$Doublet
#labels$cell_sample = data_harmony_run_label$cell_sample

labels$sample_cell = data_harmony_run_label$sample_cell

str = ''
#write.table(labels, file = paste0(path,'Labels',str,'.tsv'), quote=FALSE, sep='\t', row.names= F)
write.csv(labels, file = paste0(path,'Labels',str,'.csv'), quote=FALSE, row.names= F)


tmp1 = colnames(data_matrix)
tmp2 = labels$Cell

# metadata

sample = as.character(data_harmony_run$sample)
ident = as.character(Idents(data_harmony_run))
UmapCoord = data_harmony_run@reductions[["umap"]]@cell.embeddings
cell =  colnames(data_harmony_run)
label = as.character(Idents(data_harmony_run_label))

output = cbind(cell,sample,ident,label,UmapCoord)

dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/data.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)


# CellPhoneDB
output = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/CellPhoneDB/out baseline/'
setwd(output)
dot_plot(selected_rows = NULL,
         selected_columns = c('CD8+ T Cell|IFN+ Mono'),
         filename = 'plot.pdf',
         width = 8,
         height = 10,
         means_path = paste0(output,'/means.txt'),
         pvalues_path = paste0(output,'/pvalues.txt'),
         means_separator = '\t',
         pvalues_separator = '\t',
         output_extension = '.pdf'
)


## Seurat DE

ident1 = '8_11'
ident2 = '0'
ident1 = 'ILC1'
ident2 = 'TEMRA'
sort(unique(as.character(Idents(data_harmony_run_label))))

Features = FindMarkers(data_harmony_run_label, ident.1 = ident1, ident.2 = ident2
                       ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)

Features = Features[Features$p_val_adj < 0.05,]

path = paste0(filepath_cluster,'DE/')
dir.create( path, recursive = TRUE)
path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
print(path)
write.csv(Features, file = path,row.names=TRUE)

## Plot only clusters with gene Above threshold

cluster_list = unique(Idents(data_harmony_run_label))
cluster = cluster_list[1]

gene = 'KLRC1'

cluster_high = c()
data = data_harmony_run_label@assays[["RNA"]]@data
for (cluster in cluster_list){
  print(cluster)
  data_subset = data[,Idents(data_harmony_run_label) == cluster]
  
  if (any(data_subset[gene,] > 4) ){
    cluster_high = c(cluster_high,cluster)
  }
}

data_subset = data_harmony_run_label[,Idents(data_harmony_run_label) %in% cluster_high]

pathName <- paste0(filepath_cluster,'Plots/Subset/','High_',gene,'.png')
png(file=pathName,width=3000, height=2000,res = 100)

plot = DimPlot(data_subset,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)
dev.off()


# Save cell ids


sample_list = unique(data_harmony_run_label$sample)

for (i in 1:length(sample_list)){
  sample = sample_list[i]
  data_subset =data_harmony_run_label[,data_harmony_run_label$sample == sample]
  barcode = as.character(colnames(data_subset))
  barcode = gsub("-.*",'', barcode)
  barcode
  base = '/disk2/Projects/EloRD/Data/Harmony_Barcodes/'
  file = paste0(base, sample,'_barcode.tsv')
  write.table(barcode,file, sep='\t', row.names = F, col.names = F, quote = F)
  
}
