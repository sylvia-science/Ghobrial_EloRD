# merge the samples

# select highly variable genes

# run PCA and then 

library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(harmony)
library(ggplot2)
library(SoupX)
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
#filename_metaData = '/disk2/Projects/EloRD_Nivo_PBMC/MetaData/metaData_EloRD_Nivo_PBMC.xlsx'
metaData_old = read_excel(filename_metaData)
metaData_old = metaData_old[metaData_old$Run== 1,]


# Individual sample parameters
#sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]



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
#EloRD data + our PBMCs + PBMCs for HCA
# sample_type = 'Harmony_AllSamples_PBMC_NPBMC_Sample_Kit'
# folder_name = 'AllSamples_PBMC'
# harmony_groupby = 'Sample_Kit'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

#EloRD data + our PBMCs + PBMCs for HCA + PBMCS/BMs from HCL
# sample_type = 'AllSamples_PBMC_NPBMC_HCL_Sample_kit_tech'
# folder_name = 'AllSamples_PBMC_NPBMC_HCL'
# harmony_groupby = 'Sample_kit_tech'
# base = '/disk2/Projects/EloRD/Output/Harmony/'
# soup_base = '/disk2/Projects/EloRD/Output/Soup_MT_C100/'


#EloRD data + Nivo + Oksana data
# folder_name = 'EloRD_Nivo_Oksana'
# harmony_groupby = 'Sample_kit'
# sample_type = 'EloRD_Nivo_Oksana_Sample_kit'
# base = '/disk2/Projects/EloRD_Nivo_Oksana/Harmony/'
# filename_metaData = '/disk2/Projects/EloRD_Nivo_Oksana/MetaData/10X Sequenced Samples.xlsx'
# soup_base = '/disk2/Projects/EloRD_Nivo_Oksana/Output/Soup_MT_C100/'

# EloRD data + Nivo + Oksana data + PBMC +NPBMC - post treatment nivo
folder_name = 'EloRD_Nivo_Oksana_PBMC'
harmony_groupby = 'Sample'
sample_type = 'EloRD_Nivo_Oksana_PBMC_Sample'
base = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC/Harmony/'
filename_metaData = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC/MetaData/MetaData.csv'
soup_base = '/disk2/Projects/EloRD_Nivo_Oksana/Output/Soup_MT_C100/'
#filename_sampleParam_combine = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC/MetaData/Parameters.csv'

# EloRD data + Nivo + Oksana data + PBMC -NPBMC - post treatment nivo

folder_name = 'EloRD_Nivo_Oksana_PBMC_noNPBMC'
harmony_groupby = 'Sample_kit'
sample_type = 'EloRD_Nivo_Oksana_PBMC_Sample_noNPBMC'
base = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/'
filename_metaData = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/MetaData/MetaData.csv'
soup_base = '/disk2/Projects/EloRD_Nivo_Oksana/Output/Soup_MT_C100/'
#filename_sampleParam_combine = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC/MetaData/Parameters.csv'


# Parameters for integrated data
filename_sampleParam_combine <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam_combine)

# Load metadata
metaData = read.csv(filename_metaData)
metaData = metaData[metaData$Run== 1,]
metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]
#metaData = metaData[metaData$Sample !='BatchF',] # BatchF is crashing at quickCluster in scran norm

# Combine new metadata with old metadata
#metaData$Sample[!(metaData$Sample %in% metaData_old$Sample)]
#metaData_old$Sample[!(metaData_old$Sample %in% metaData$Sample)]
#tmp = metaData_old[!(metaData_old$Sample %in% metaData$Sample),]

metaData_new = merge(metaData,metaData_old, by = 'Sample', all.x = TRUE,sort = FALSE)
metaData = metaData_new[match( metaData$Sample,metaData_new$Sample),]
names(metaData) = gsub(".x","",names(metaData),fixed = TRUE)


# Get params
PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 

folder = paste0(base,
                '/Batch_',harmony_groupby,'/')

#folder = paste0(base,folder_name,
#                '/Batch_',harmony_groupby,'/')
dir.create(folder,recursive = T)
run = F




sample_list = as.character(metaData$Sample)
#sample_list[!(sample_list %in% unique(data_harmony_run$sample))]

# For running individual sample
i = 1
sample_name = 'BatchF'

unique(data_harmony_run$`Patient Number`)


if (run){
  ################################
  # Integrate data with Harmony
  ################################
  
  # Put filtered and normalized data into list
  data_list = vector(mode = "list",length = length(sample_list))
  data_list_norm = vector(mode = "list",length = length(sample_list))
  #for (i in 1:nrow(sampleParam)){
  for (i in 1:length(sample_list)){
    sample_name = sample_list[i]
    print(sample_name)
    # Load from soup data filtered for MT and count
    
    folder_input = paste0('/disk2/Data/Soup_MT_Count/', sample_name , '/')
    
    #/disk2/Projects/EloRD/Output/Soup_MT_C100/
    data_i = loadRData(paste0(folder_input,sample_name,'.Robj'))
    data_i$sample = sample_name

    # remove empty drops
    data_i = load_emptyDrops(data_i, '/disk2/Data/EmptyCells/')
    data_i$is_cell[is.na(data_i$is_cell) ] = T # NA values are ones that didn't pass the count threshold
    print(sum(!data_i$is_cell == T))
    data_i = data_i[,data_i$is_cell]
    
    # Remove features below threshold  
    min_feature = 100
    data_i = data_i[,data_i$nFeature_RNA > min_feature]
  
    
    # Remove MT > 15 already done
    print(ncol(data_i))
    if (ncol(data_i) > 100){ # Scran Normalization will crash with less cells
      data_list[[i]] = data_i
      data_list_norm[[i]] = ScranNorm(data_i)
      
    }
    
  }
  
  data_list_norm =data_list_norm[lengths(data_list_norm) != 0]
  data_merge = merge(x =  data_list_norm[[1]] ,y = data_list_norm[2:length(data_list_norm)], merge.data = T)
  
  data_merge = addMetaData(data_merge, metaData)
  #data_merge = load_emptyDrops(data_merge, '/disk2/Data/EmptyCells/')
  data_merge = load_Doublets(data_merge)
  data_merge$kit = data_merge$`10X kit`
  data_merge$split_var = ''
  cell_names_all = sub("_.*", "", colnames(data_merge)) 
  # Create new unique cell names so we can match them to individual samples
  data_merge$cell_sample = paste0(data_merge$sample ,' ',cell_names_all) 
  data_merge$sample_cell = paste0(cell_names_all,' ',data_merge$sample )
   
  data_merge_run = FindVariableFeatures(data_merge, selection.method = "vst", nfeatures = 2000)
  #data_harmony_run = ScaleData(data_harmony_run, vars.to.regress = c('kit'))
  data_merge_run = ScaleData(data_merge_run)
  
  # Run PCA and visualize
  PCA_dim_init = 90
  data_merge_run = RunPCA(data_merge_run,npcs = PCA_dim_init)
  reduction = 'pca'
  visualize_PCA(data_merge_run,folder,PCA_dim_init,reduction)
  
  
  #PCA_dim = 80
  #resolution_val = 1
  # Run pipeline on merged data
  data_merge_run = RunUMAP(reduction = "pca",data_merge_run, dims = 1:PCA_dim)
  data_merge_run = FindNeighbors(data_merge_run, reduction = "pca", dims = 1:PCA_dim)
  data_merge_run = FindClusters(data_merge_run,resolution = resolution_val)
  
  data_merge_run$FeatureLessThan400 = data_merge_run$nFeature_RNA < 400
  data_merge_run$kit = data_merge_run$X10X.kit
  
  # Save merged data
  folder = paste0(base,
                  '/Merge/')
  dir.create(folder)
  path = paste0(folder,'data_merge_run','.Robj')
  save(data_merge_run,file= path)
  #data_merge_run = loadRData((path))
  
  data_merge_run$split_var = ''
  
  # AddMetaData again if it has changed from last run
  #data_merge_run = addMetaData(data_merge_run, metaData)
  #data_merge_run$kit = data_merge_run$`10X kit`
  # Plot merged data
  cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
  groupBy_list = c('sample','Diagnosis','kit',
                   'Batch','LowCount',
                   'FeatureLessThan400','Sample Type','Gender','Race')
  featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
  splitBy_list = NA
  
  # DR calc will crash with too many cells
  #data = as.data.frame(data_merge_run@assays[["RNA"]]@counts)
  #DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
  #data_merge_run$DR = DR
  

  plotAll(data_merge_run, folder = folder,
          sample_name,sampleParam = NA,
          cell_features = cell_features, plot_PCA = F,
          label_TF = F,  DE_perm_TF = F, 
          markersTF = F,
          groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
          PCA_dim = PCA_dim,resolution_val = resolution_val)
  
  
  # Run harmony and pipeline
  harmony_dim = PCA_dim
  data_harmony_run = RunHarmony(data_merge_run,group.by.vars =  c('sample','kit'),
                                dims.use = 1:harmony_dim)
  
  folder = paste0(base,'/Batch_',harmony_groupby,'/')
  reduction = 'harmony'
  visualize_PCA(data_harmony_run,folder,harmony_dim,reduction)
  
  data_harmony_run = RunUMAP(reduction = "harmony",data_harmony_run, dims = 1:harmony_dim)
  data_harmony_run = FindNeighbors(data_harmony_run, reduction = "harmony", dims = 1:harmony_dim)
  data_harmony_run = FindClusters(data_harmony_run,resolution = resolution_val)
  
  
  plot = DimPlot(data_harmony_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
  # Save integration
  path = paste0(folder,'data_run','.Robj')
  save(data_harmony_run,file= path)
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  dir.create(filepath_cluster,recursive = T)

}else{
  #######################
  # Load processed data
  #######################
  #resolution_val = 3
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )

  # Load data and change resolution
  path = paste0(folder,'data_run','.Robj')
  data_harmony_run = loadRData(path)
  tmp = data_harmony_run@meta.data[paste0('RNA_snn_res.', resolution_val)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_harmony_run) = tmp
  data_harmony_run$orig.ident = Idents(data_harmony_run)
  
  # Reload vars in case they have changed
  cell_names_all = sub("_.*", "", colnames(data_harmony_run))
  data_harmony_run$cell_sample = paste0(cell_names_all,' ',data_harmony_run$sample )
  data_harmony_run$sample_cell = paste0(data_harmony_run$sample,' ',cell_names_all )
  #data_harmony_run = load_Doublets(data_harmony_run)
  data_harmony_run$FeatureLessThan400 = data_harmony_run$nFeature_RNA < 400
  data_harmony_run = addMetaData(data_harmony_run, metaData)
  data_harmony_run$kit = data_harmony_run$X10X.kit
  
  
  #cell_names_main = as.character(data_harmony_run$cell_sample)

  path = '/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Cluster/PCA30/res3//Data/Labels.csv'
  data_harmony_run = addOldLabels(path,data_harmony_run, 'OldCellType')

  data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)


  if (sample_type == 'Harmony_AllSamples_Sample_Kit'){
    data_harmony_run$sample_cell = paste0(data_harmony_run$Sample,'_', colnames(data_harmony_run))
  
    TSCM_res4 = read.csv(paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Subcluster/T Cell/Cluster/PCA30/res3.5/'
                                ,'Data/','TSCM_res4','.csv'))
  
    new_Idents = as.character(Idents(data_harmony_run_label))
    new_Idents[colnames(data_harmony_run_label) %in% TSCM_res4$x] = 'TSCM'
    Idents(data_harmony_run_label) = new_Idents
  }
  
}

## Add subset labels
#######################################
## Original EloRD run

Ident_list = as.character(Idents(data_harmony_run_label))
path = '/home/sujwary/Desktop/scRNA/Output/Harmony//Batch_Sample_Kit/Subcluster/T Cell/Cluster/PCA30/res3.5//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'TcellLabels')
Ident_list[data_harmony_run_label$TcellLabels!=''] = as.character(data_harmony_run_label$TcellLabels[data_harmony_run_label$TcellLabels!=''])

        
path = '/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Subcluster/Mono_DC/Cluster/PCA30/res1.6/Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'monoDCLabels')
Ident_list[data_harmony_run_label$monoDCLabels!=''] = as.character(data_harmony_run_label$monoDCLabels[data_harmony_run_label$monoDCLabels!=''])

path = '/home/sujwary/Desktop/scRNA/Output/Harmony//Batch_Sample_Kit/Subcluster/NK_RemoveRiboFeatures/Cluster/PCA30/res3//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'NKLabels')
Ident_list[data_harmony_run_label$NKLabels!=''] = as.character(data_harmony_run_label$NKLabels[data_harmony_run_label$NKLabels!=''])
Idents(data_harmony_run_label) = Ident_list

## EloRD_Nivo_Oksana
Ident_list = as.character(Idents(data_harmony_run_label))
path = '/disk2/Projects/EloRD_Nivo_Oksana/Harmony/Batch_Sample_kit/Subcluster/T Cell/Cluster/PCA30/res3/Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'TcellLabels')
Ident_list[data_harmony_run_label$TcellLabels!=''] = as.character(data_harmony_run_label$TcellLabels[data_harmony_run_label$TcellLabels!=''])

path = '/disk2/Projects/EloRD_Nivo_Oksana/Harmony//Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2/Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'monoDCLabels')
Ident_list[data_harmony_run_label$monoDCLabels!=''] = as.character(data_harmony_run_label$monoDCLabels[data_harmony_run_label$monoDCLabels!=''])

path = '/disk2/Projects/EloRD_Nivo_Oksana/Harmony//Batch_Sample_kit/Subcluster/NK/Cluster/PCA30/res3.5//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'NKLabels')
Ident_list[data_harmony_run_label$NKLabels!=''] = as.character(data_harmony_run_label$NKLabels[data_harmony_run_label$NKLabels!=''])
Idents(data_harmony_run_label) = Ident_list

## EloRD_Nivo_Oksana_PBMC_noNPBMC
Ident_list = as.character(Idents(data_harmony_run_label))
path = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony//Batch_Sample_kit/Subcluster/T Cell/Cluster/PCA30/res3.5//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'TcellLabels')
Ident_list[data_harmony_run_label$TcellLabels!=''] = as.character(data_harmony_run_label$TcellLabels[data_harmony_run_label$TcellLabels!=''])

path = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony//Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'monoDCLabels')
Ident_list[data_harmony_run_label$monoDCLabels!=''] = as.character(data_harmony_run_label$monoDCLabels[data_harmony_run_label$monoDCLabels!=''])

path = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony//Batch_Sample_kit/Subcluster/NK/Cluster/PCA30/res2//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'NKLabels')
Ident_list[data_harmony_run_label$NKLabels!=''] = as.character(data_harmony_run_label$NKLabels[data_harmony_run_label$NKLabels!=''])

path = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony//Batch_Sample_kit/Subcluster/Stem_Cell_B_Cell_RegressCellCycle/Cluster/PCA30/res3//Data/labels.csv'
data_harmony_run_label = addOldLabels(path,data_harmony_run_label, 'HSCLabels')
Ident_list[data_harmony_run_label$HSCLabels!=''] = as.character(data_harmony_run_label$HSCLabels[data_harmony_run_label$HSCLabels!=''])
Idents(data_harmony_run_label) = Ident_list

#############################################
#metaData = read_excel(filename_metaData)


#table(data_harmony_run$sample)
#table(data_harmony_run$Treatment)
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
dir.create(filepath_cluster,recursive = T)


##
#data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)
##

##########################
## Clean up labeled data
##########################
remove_list = c('Remove','14','32','Erythrocyte','DC/T-Cell DBL',
                'Mono/CD8+ T Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL','42','Plasma Cell',0:60,
                'CD14+ Mono/T-cell DBL','CD14+ Mono/CD8+ T-cell DBL','CD16+ Mono/T-cell DBL','dMono','dNK','Plasma cell','Pro Erythrocyte',
                'DC/T-cell DBL','MK')

remove_list = c('Remove','Erythrocyte','DC/T-Cell DBL',
                'Mono/CD8+ T Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL',1:60,
                'CD14+ Mono/T-cell DBL','CD14+ Mono/CD8+ T-cell DBL','CD16+ Mono/T-cell DBL',
                'dMono','dNK','Pro Erythrocyte','Eryrhrocyte',
                'DC/T-cell DBL','MK')

remove_list = c('Remove','14','32','Erythrocyte',
                '42','Plasma Cell','41','',0:100,
                'Remove 29','Remove 30')

remove_list = c(0:100,'Remove')

data_harmony_run_label_remove = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
Idents(data_harmony_run_label_remove) = as.character(Idents(data_harmony_run_label_remove))
sort(unique(Idents(data_harmony_run_label_remove)))
#data_harmony_run_label_remove = data_harmony_run_label[,data_harmony_run_label$`Sample Type` == 'PBMC']


#data_harmony_run_label = data_harmony_run_label[,data_harmony_run_label$is_cell]
#data_harmony_run_label = data_harmony_run_label[,Idents(data_harmony_run_label) != 'Remove']


plot = DimPlot(data_harmony_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 3)
print(plot)

plot = DimPlot(data_harmony_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 2)
print(plot)

plot = DimPlot(data_harmony_run_label_remove,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

# Plot old labels
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

groupBy_list = c('sample','Diagnosis','kit','Gender',
                 'Treatment','Batch','LowCount',
                 'Doublet','FeatureLessThan400',
                 'Sample Type','Race')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA

#data = as.data.frame(data_harmony_run@assays[["RNA"]]@counts)
#DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
#data_harmony_run$DR = DR
unique(data_harmony_run$`Sample Type`)
unique(data_harmony_run$sample[data_harmony_run$`Sample Type` == 'Peripheral Blood'])

plotAll(data_harmony_run, folder = folder, 
        sample_name,sampleParam = NA,
        cell_features = cell_features, plot_PCA = F,
        label_TF = F,  DE_perm_TF = F, DE_perm_clusters = NA,
        markersTF = T,
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val)


plotAll(data_harmony_run_label, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,plot_PCA = F,
        label_TF = F, DE_perm_TF = F, DE_perm_clusters = NA,
        markersTF = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val,str = '_label')

data_harmony_run_label_remove

plotAll(data_harmony_run_label_remove, folder = folder,
        sample_name,sampleParam = NA,
        cell_features = cell_features,plot_PCA = F,
        label_TF = F,  DE_perm_TF = F, DE_perm_clusters = NA,
        markersTF = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim,resolution_val = resolution_val,str = '_label_remove')

###################################
## Get old cell type per cluster
###################################
dir.create(paste0(filepath_cluster,'/Data'))
cell_table = table( data_harmony_run$OldCellType,Idents(data_harmony_run))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','OldCellTypeByIdent_Table.csv'))

cell_table = table( data_harmony_run_label$sample,Idents(data_harmony_run_label))
cell_table_new = as.data.frame.matrix(cell_table)
cell_table_new$Sample = rownames(cell_table_new)
#cell_table_new = merge(cell_table_new,metaData[,c('Sample','Diagnosis','Study','Treatment')], by = 'Sample', all.x = TRUE,sort = FALSE)
#cell_table_new = merge(cell_table_new,metaData_old[,c('Sample','Treatment')], by = 'Sample', all.x = TRUE,sort = FALSE)
#cell_table_new = cell_table_new[,!(colnames(cell_table_new) %in% remove_list)]
write.csv(cell_table_new, file = paste0(filepath_cluster,'/Data/','sampleByIdent_Table.csv'))

write.csv(metaData, file = paste0(filepath_cluster,'/Data/','metadata.csv'))


sort(colnames(cell_table_new))

cell_table = table( as.character(data_harmony_run$sample))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','NumCellPerSample.csv'))

cell_table = table( data_harmony_run_label$sample,Idents(data_harmony_run_label))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','ClusterPerSample.csv'))
tmp = table(metaData$Diagnosis,metaData$Study )


cell_table = table( data_harmony_run_label_remove$sample,Idents(data_harmony_run_label_remove))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','ClusterPerSample_subcluster_label.csv'))

Idents = as.character(Idents(data_harmony_run_label))
Idents[Idents %in% 1:100] = 'Remove'
Idents(data_harmony_run_label) =  Idents
plot = DimPlot(data_harmony_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 2)
print(plot)
cell_table = table( data_harmony_run_label$sample,Idents(data_harmony_run_label))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','ClusterPerSample_all.csv'))

metadata_sample = unique(metaData$Sample)
data_sample = unique(data_harmony_run$sample)
metadata_sample[!(metadata_sample %in% data_sample)]
data_sample[!(data_sample %in% metadata_sample)]


for (cluster in sort(unique(Idents(data_harmony_run)))){
  print(cluster)
  data_harmony_run_subset = data_harmony_run[,Idents(data_harmony_run) == cluster]
  print(table(data_harmony_run_subset$OldCellType))
}


########################
## Plot markers
########################
cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')
#cell_features = cell_features[cell_features$Cell == 'monocyte_FCGR3A',]
filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )

PlotKnownMarkers(data_harmony_run, folder = paste0(filepath_cluster,'Cell Type/FeaturePlotFix/'), cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F, plotAll = F)

cell_features = cell_features[cell_features$Plot_marker == 1,]
PlotKnownMarkers(data_harmony_run, folder = paste0(filepath_cluster,'Cell Type/AllGenes/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F, plotAll = T) # Plot all genes without cell type prefixes
print(DimPlot(data_harmony_run, reduction = "umap", label = TRUE, pt.size = .1))

gene_list = c('IRF4','IRS1','CARM1','CFLAR','PRDM1','MARCH5','FURIN','HERPUD1')

gene_list = c('CRBN')
FeaturePlot_GeneList(data_harmony_run,gene_list,
                     folder = paste0(filepath_cluster,'Cell Type/Mahshid/'),
                     FeaturePlotFix = T,str = '')

StackedVlnPlotHelper(data_harmony_run_label,gene_list,
                                folder_heatMap = paste0(filepath_cluster,'Plots/Violin/'),
                                filename = paste0('CRBN_AllCluster.pdf'),width = 20)
 
# Plot violin plot per cluster
celltype_list = unique(Idents(data_harmony_run_label_remove))
pdf(paste0(filepath_cluster,'Plots/Violin/CRBN.pdf'))
for (celltype in celltype_list){
  print(celltype)
  data_subset = data_harmony_run_label_remove[,Idents(data_harmony_run_label_remove) == celltype]
  plot= VlnPlot(data_subset, 'CRBN', pt.size = 0.1) 
  print(plot)
}
dev.off()
##################################
## Violin plot for specific gene
##################################
path = '/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Subcluster/Mono_DC/Cluster/PCA30/res1.6/Data/labels.csv'
celltype = 'MonoDC'
path = '/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Subcluster/NK_RemoveRiboFeatures/Cluster/PCA30/res3/Data/labels.csv'
celltype = 'NK'
path  = '/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Subcluster/T Cell/Cluster/PCA30/res3.5/Data/labels.csv'
celltype = 'TCell'

label = read.csv(path)
data_subset_gene = t(data_harmony_run_label_remove['CRBN',][["RNA"]]@data)
data_subset_gene = as.data.frame(data_subset_gene)
data_subset_gene$label = Idents(data_harmony_run_label_remove)
data_subset_gene$sample = data_harmony_run_label_remove$sample
data_subset_gene$sample_cell = data_harmony_run_label_remove$sample_cell
data_subset_gene$Treatment =  data_harmony_run_label_remove$Treatment                 
                                               
data_subset_gene = data_subset_gene[data_subset_gene$sample_cell %in% label$sample_cell,]

pdf(paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Cluster/PCA30/res3/Plots/Violin/','CRBN_',celltype,'.pdf'), width =12)
plot = ggplot(data_subset_gene) + geom_violin(aes(x=label,y=CRBN,fill=Treatment)) + theme_bw() +
  theme(axis.text=element_text(size=12),  axis.title=element_text(size=14,face="bold")) + 
  theme(axis.text=element_text(angle=65,hjust=1))
print(plot)                                               
dev.off()


unique(data_subset_gene$label)
write.csv(data_subset_gene, 
      paste0('/home/sujwary/Desktop/scRNA/Output/Harmony/Batch_Sample_Kit/Cluster/PCA30/res3/Data/',celltype,'_CRBN.csv'))


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



######################################
### Plot DE genes in each cluster
######################################
str = ''
str = '_label'
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

##########################
### Plot DE permutations
##########################
data_input = data_run_subset_label
data_input = data_harmony_run
cluster_list1 = sort(unique(Idents(data_input)))
cluster_list2 = cluster_list1

cluster_list1 = c('CMPC','HSC','MDPC', 'GMPC','Pro B-cell', 'Pre B-cell')
cluster_list2 = c('CMPC','HSC','MDPC', 'GMPC','Pro B-cell', 'Pre B-cell')


cluster1 = 27
cluster2 = 24
str = '_label_remove'
str = '_label'
#str = ''
for (cluster1 in cluster_list1){
  #
  cluster1 = gsub("/", "", cluster1, fixed = T)
  filepath_mvn = paste0( filepath_cluster, 'DE/nVsm/')
  filepath = paste0(filepath_mvn
                    ,'Features_',cluster1,'Vsn',str,
                    '.csv')
  DE_cluster = read.csv(file = filepath)
  
  for (cluster2 in cluster_list2){
    #data_cluster  = data_run_subset[, Idents(data_run_subset) == cluster]
    #print(cluster)
    #print(summary(data_cluster$origIdent))
    cluster2 = gsub("/", "", cluster2, fixed = T)
    DE_cluster_m = DE_cluster[DE_cluster$ident_2 == cluster2,]
    DE_cluster_m = DE_cluster_m[order(-abs(DE_cluster_m$avg_logFC )),]
    DE_cluster_m = DE_cluster_m[DE_cluster_m$p_val_adj < 0.05 & DE_cluster_m$avg_logFC > 0 ,]
    DE_cluster_m = DE_cluster_m[!grepl("MT-", DE_cluster_m$gene),]
    DE_cluster_m = DE_cluster_m[!grepl("^RP[SL]", DE_cluster_m$gene),]
    
    gene_list = DE_cluster_m[1:20,]
    gene_list = (gene_list[!rowSums(is.na(gene_list)) == ncol(gene_list),])
    gene_list$Markers = gene_list$gene
    #colnames(gene_list)[grepl("gene", colnames(gene_list))] <- "Markers"
    
    #next
    if (nrow(gene_list) > 0){
      gene_list$Cell = paste0(cluster1,'Vs',cluster2)
      gene_list$Plot_marker = 1
      path = paste0(filepath_mvn,'FeaturePlot/Cluster',cluster1,str,'/')
      PlotKnownMarkers(data_input, 
                       folder = path, 
                       cell_features = gene_list,
                       plotType ='FeaturePlotFix' , 
                       str = '',markerSize = 2, prefix_logFC = T)
    }
  }
  
  
}

##################################
## Heatmap for specific genes
##################################


data_harmony_run_label = ScaleData(data_harmony_run_label,features =  rownames(data_harmony_run_label))

remove_list= c('HSC','Pro Erythrocyte','MDPC','MK','Erythrocyte','Pro B Cell','Pre B Cell','Remove',
  'CMPC','42','GMPC','Plasma Cell','B Cell')

data_harmony_run_label_remove = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]

gene_list = gene_list = c('AXL','GAS6')
folder = paste0(filepath_cluster,'/Plots/Heatmap/')
dir.create(folder,recursive = T)
pathName <- paste0(folder,paste0('HeatMap', '_AXL_GAS6','.pdf'))
pdf(file=pathName, height=12, width=20)
#png(file=pathName)
#png(file=pathName,width=2000, height=1000,res = 100)

plot = DoHeatmap(data_harmony_run_label_remove, features = gene_list)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()

########################
## Entropy
#########################
cluster_list = levels(unique(Idents(data_harmony_run)))

for (cluster in cluster_list){
  data_subset = data_harmony_run[,Idents(data_harmony_run) == cluster]
  print(cluster)
  cell_percent =  100*table(data_subset$GeneralCellType)/ncol(data_subset)
  print(cell_percent)
  print(max(cell_percent))
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
## Empty Cells
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
#######################################
## NMF with specific cell types
#######################################
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

celltype = 'all'


unique(Idents(data_harmony_run_label_remove))


cell_list[!(cell_list %in% unique(Idents(data_harmony_run_label_remove)))]

unique(Idents(data_harmony_run_label_remove))[!(unique(Idents(data_harmony_run_label_remove)) %in%cell_list )]

#data_input = data_harmony_run_label_remove[,Idents(data_harmony_run_label_remove) %in% cell_list]
data_input = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in%
                                         c( 0:100,'Remove', 'B-cell','Pro B-cell','Pre B-cell','Erythrocyte','Pro Erythrocyte',
                                            'Neutrophil','MK','HSC','CMPC','MDPC', 'Plasma Cell','pDC','cDC1','cDC2'))]

data_input = data_harmony_run_label

plot = DimPlot(data_input,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_matrix = data_input@assays[["RNA"]]@data
data_matrix_var = data_matrix[rownames(data_matrix) %in% data_input@assays[["RNA"]]@var.features,]

path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_PreRerun','.tsv')
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')

path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_varFeature_subset_PreRerun','.tsv')
write.table(data_matrix_var, 
            file=path, 
            quote=FALSE, sep='\t')


# UmapCoord and metadata

#data_input = data_harmony_run_label
sample = as.character(data_input$sample)
ident = as.character(Idents(data_harmony_run[,colnames(data_harmony_run_label)]))
UmapCoord = as.data.frame(data_input@reductions[["umap"]]@cell.embeddings)
cell =  colnames(data_input)
label = as.character(Idents(data_input))
old_ident = as.character(data_input$OldCellType)
#old_ident_T = as.character(data_input$OldTcellLabels)
#old_ident_NK = as.character(data_input$OldNKLabels)
#old_ident_Mono = as.character(data_input$OldMonoLabels)


output = cbind(cell,sample,label,UmapCoord,ident)
output$Patient = data_input$Patient
output$Diagnosis = data_input$Diagnosis
output$Study = data_input$Study
output$Treatment = data_input$Treatment
output$Response = data_input$Response 
#output = output[,c('sample')]
#output = data_run_subset_label@meta.data
dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/data_',celltype,'_PreRerun.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)

## Rerun data and Save

harmony_dim = 60
data_input_run = FindVariableFeatures(data_input, selection.method = "vst", nfeatures = 2000)
data_input_run = RunUMAP(reduction = "harmony",data_input_run, dims = 1:harmony_dim)
data_input_run = FindNeighbors(data_input_run, reduction = "harmony", dims = 1:harmony_dim)
data_input_run = FindClusters(data_input_run,resolution = 3)

data_input_run$ident_orig = Idents(data_input)
plot = DimPlot(data_input_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

plot = DimPlot(data_input_run,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8,group.by  = 'ident_orig')
print(plot)

path = paste0(folder,'data_',celltype,'_run','.Robj')
save(data_input_run,file= path)
####################
sample = as.character(data_input_run$sample)
ident = as.character(Idents(data_input_run))
UmapCoord = data_input_run@reductions[["umap"]]@cell.embeddings
ident_orig = Idents(data_input)

output = as.data.frame(sample)
output$ident = ident
output$ident_orig = data_input_run$ident_orig

output = cbind(output,UmapCoord)
#output = data_harmony_run_label@meta.data

dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/data_',celltype,'_Rerun','.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)

data_matrix = data_input_run@assays[["RNA"]]@data
data_matrix_var = data_matrix[rownames(data_matrix) %in% data_input_run@assays[["RNA"]]@var.features,]
path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_varFeature_Rerun','.tsv')
write.table(data_matrix_var, 
            file=path, 
            quote=FALSE, sep='\t')

################################
# Save all counts in matrix for nmf
################################
path = paste0('/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit',
              '.tsv')
data_matrix = data_harmony_run@assays[["RNA"]]@data
data_matrix = data_matrix[rownames(data_matrix) %in% data_harmony_run@assays[["RNA"]]@var.features,]
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')


#data_harmony_run_label = data_harmony_run_label[,Idents(data_harmony_run_label) %in% c('NK','CD14+ Mono')]
#data_harmony_run_label = label_cells(data_harmony_run,cluster_IDs)

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



########################
# Save umap and metadata
#########################
sample = as.character(data_harmony_run$sample)
ident = as.character(Idents(data_harmony_run))
UmapCoord = data_harmony_run@reductions[["umap"]]@cell.embeddings
cell =  colnames(data_harmony_run)
label = as.character(Idents(data_harmony_run_label))

output = cbind(cell,sample,ident,label,UmapCoord)

dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/data.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)

#########################
# CellPhoneDB
#########################
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

##############################
## Seurat DE
##############################
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



