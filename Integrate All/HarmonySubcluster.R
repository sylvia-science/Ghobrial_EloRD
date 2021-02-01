# merge the samples

# select highly variable genes

# run PCA 

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
library(edgeR)
library(SingleCellExperiment)
library(MAST)
library(UpSetR)
library(SingleCellSignalR)
library(fgsea)
#data(examplePathways)
#data(exampleRanks)



source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/FunctionsIntegrate.R')
source('/home/sujwary/Desktop/scRNA/Code/LoadCellData.R')
source('/home/sujwary/Desktop/scRNA/Code/Integration/PlotAll.R')
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Integrate All/Entropy.R')
source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')

celltype = 'T Cell'
cell_list = c('T-cell','CD8+ T-cell')
#resolution_val_subset = 3.5
# # 
celltype = 'Mono_DC'
cell_list = c('CD14+ Mono','CD16+ Mono','DC')
resolution_val_subset = 1.6
# # # #
# # celltype = 'NK_Remove_Nfeature2000'
# # cell_list = c('NK')
# # resolution_val_subset = 1.2
# # 
celltype = 'NK_RemoveRiboFeatures'
cell_list = c('NK')
resolution_val_subset = 3


# 
# celltype = 'NK'
# cell_list = c('NK')
# resolution_val_subset = 3
# # # 
# celltype = 'B_Cell'
# cell_list = c('B-cell','Pro B-cell','Pre B-cell','Plasma Cell')
# resolution_val_subset = 3

# celltype = 'Stem_Cell'
# cell_list = c('HSC','CMPC', 'GMPC','MDPC','42')
# resolution_val_subset = 3


# celltype = 'Remove_cluster'
# cell_list = c('Remove','9')
# resolution_val_subset = 3


filename_sampleParam = paste0('/home/sujwary/Desktop/scRNA/Param/','sample_parameters_Scran.xlsx')
sampleParam <- read_excel(filename_sampleParam)
filename = paste0('/home/sujwary/Desktop/scRNA/Param/','Cluster_ID_testNorm.xlsx')
cluster_id_param = read_excel(filename)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]
#metaData = metaData[metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]


sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

#downsample = read.csv('/home/sujwary/Desktop/scRNA/Output/CompareIntegration/DownSampleCells.csv')
downsample  = NA #downsample$x


sample_type = 'Harmony_AllSamples_Sample_Kit'
#sample_type = 'Harmony_PBMC_Sample_Kit'
folder_name = 'AllSamples'
base = '/home/sujwary/Desktop/scRNA/Output/Harmony/'
harmony_groupby = '_Sample_Kit'
# folder_name = 'AllSamples'
# #folder_name = 'PBMC'
# base = '/home/sujwary/Desktop/scRNA/Output/Harmony/'
# 
# sample_type = 'Harmony_AllSamples_PBMC_NPBMC_Sample_Kit'
# folder_name = 'AllSamples_PBMC'
# harmony_groupby = 'Sample_Kit'
# base = '/disk2/Projects/EloRD/Output/Harmony/'
# 
# sample_type = 'AllSamples_PBMC_NPBMC_HCL_Sample_kit_tech'
# folder_name = 'AllSamples_PBMC_NPBMC_HCL'
# harmony_groupby = '_Sample_kit_tech'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

# sample_type = 'AllSamples_PBMC_NPBMC_HCL_Sample_kit'
# folder_name = 'AllSamples_PBMC_NPBMC_HCL'
# harmony_groupby = '_Sample_kit'
# base = '/disk2/Projects/EloRD/Output/Harmony/'

PCA_dim = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == sample_type]
resolution_val = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == sample_type]
cluster_IDs = sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == sample_type] 
cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
resolution_val_subset = sampleParam_combine$resolution_val[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
PCA_dim_subset = sampleParam_combine$PCA_dim[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]

harmony_dim_subset = sampleParam_combine$harmony_dim[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]


#resolution_val_subset = 1.2
#PCA_dim_subset = 20
#harmony_dim_subset = 20

#patient_list = c(12, 16, 20)

i = 1

filename_testIntRun= '/home/sujwary/Desktop/scRNA/Param/TestIntegrationRuns.xlsx'
Samples_runs = read_excel(filename_testIntRun)

#folder = 'Intra-v3_1'
#folder = 'Inter-version'

sample_list = metaData$Sample

folder_main = paste0(base,folder_name,
                '/Batch',harmony_groupby,'/')
dir.create(folder_main,recursive = T)
folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val_subset,'/' )
dir.create(filepath_cluster,recursive = T)



run = F

folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
cell_features_file = paste0('/home/sujwary/Desktop/scRNA','/Data/Cell_IDS.xlsx')
cell_features = read_excel(cell_features_file)
cell_features = cell_features[cell_features$Plot_marker == 1,]

filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',PCA_dim_subset,'/res',resolution_val_subset,'/' )

if (run){
  path = paste0(folder_main,'data_run','.Robj')
  data_orig = loadRData(path)
  tmp = data_orig@meta.data[paste0('RNA_snn_res.', resolution_val)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_orig) = tmp
  
  data_orig =label_cells(data_orig,cluster_IDs)
  data_orig = data_orig[,Idents(data_orig) %in% cell_list]

  path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK/Cluster/PCA30/res3//Data/cellsIdents.csv'
  old_cells = read.csv(path)
  remove_list = c(0, 8, 11, 17, 18, 19, 21, 22)
  #remove_list = c( 8, 11, 17, 18, 19, 21, 22)
  remove_list = c()
  
  remove_cells = old_cells[(old_cells$Ident %in% remove_list),]
  data_orig = data_orig[, !(as.character(data_orig$cell_sample)  %in% remove_cells$cell_sample)]
  
  plot = DimPlot(data_orig,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
  cell_names_all = colnames(data_orig)
  cell_names_all = sub("_.*", "", cell_names_all)
  data_orig$cell_sample = paste0(cell_names_all,' ',data_orig$sample )
  
  
  data_merge_subset =  data_orig
  
  data_merge_subset = addMetaData(data_merge_subset, metaData)
  data_merge_subset = load_emptyDrops(data_merge_subset)
  data_merge_subset = load_Doublets(data_merge_subset)
  #data_merge_subset = load_CellLabel(data_merge_subset)
  data_merge_subset$origIdent = Idents(data_orig)
  data_merge_subset$kit = data_merge_subset$`10X kit`
  data_merge_subset$split_var = ''
  
  data_run_subset = FindVariableFeatures(data_merge_subset, selection.method = "vst", nfeatures = 2000)
  
  
  # Remove Ribo genes form var genes
  var_features = data_run_subset@assays[["RNA"]]@var.features
  sum(grepl("^RP[SL]",var_features))
  sum(grepl("^MT",var_features))
  var_features = var_features[!grepl("^RP[SL]",var_features)]
  data_run_subset@assays[["RNA"]]@var.features  = var_features
  
  data_run_subset = ScaleData(data_run_subset)
  PCA_dim_orig = 50
  data_run_subset = RunPCA(data_run_subset,npcs = PCA_dim_orig)
  
  reduction = 'pca'
  visualize_PCA(data_run_subset,folder_subcluster,PCA_dim_orig,reduction)

  #PCA_dim = 20
  data_run_subset = RunUMAP(reduction = "pca",data_run_subset, dims = 1:PCA_dim_subset)
  data_run_subset = FindNeighbors(data_run_subset, reduction = "pca", dims = 1:PCA_dim_subset)
  data_run_subset = FindClusters(data_run_subset,resolution = 1)
  
  path = paste0(folder_subcluster,'data_merge_run','.Robj')
  save(data_run_subset,file= path)
  
  data_run_subset$split_var = ''
  groupBy_list = c('sample','Diagnosis','kit',
                   'Treatment','Batch','LowCount',
                  'Doublet',  'FeatureLessThan400','Technology')
  featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
  splitBy_list = NA
  
  data = as.data.frame(data_run_subset@assays[["RNA"]]@counts)
  DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
  data_run_subset$DR = DR
  
  plotAll(data_run_subset, folder = paste0(folder_subcluster,'Merge/'),
          sample_name,sampleParam = NA,
          cell_features = cell_features, plot_PCA = F,
          label_TF = F,integrate_TF = F,  DE_perm_TF = F, 
          clusterTF =F, markersTF = F,
          groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
          PCA_dim = PCA_dim_subset,resolution_val = 1)
  
  data_run_subset = RunHarmony(data_run_subset,group.by.vars =  c("sample", "10X kit",'Technology'),
                               dims.use = 1:harmony_dim_subset)
  
  #data_run_subset = RunHarmony(data_run_subset,group.by.vars =  c("sample", "10X kit"),
  #                            dims.use = 1:harmony_dim_subset)
  
  reduction = 'harmony'
  visualize_PCA(data_run_subset,folder_subcluster,harmony_dim_subset,reduction)

  
  data_run_subset = RunUMAP(reduction = "harmony",data_run_subset, dims = 1:harmony_dim_subset)
  data_run_subset = FindNeighbors(data_run_subset, reduction = "harmony", dims = 1:harmony_dim_subset)
  data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  data_run_subset$origIdent = Idents(data_orig)
  
  
  dir.create(folder_subcluster,recursive = T)
  path = paste0(folder_subcluster,'data_run_PC',harmony_dim_subset,'.Robj')
  save(data_run_subset,file= path)
  ###########################
}else{
  
  folder_subcluster = paste0(folder_main, 'Subcluster/',celltype,'/')
  filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',PCA_dim_subset,'/res',resolution_val_subset,'/' )
  
  path = paste0(folder_subcluster,'data_run_PC',harmony_dim_subset,'.Robj')
  #path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  
  cell_names_all = sub("_.*", "", colnames(data_run_subset))
  data_run_subset$cell_name = cell_names_all
  data_run_subset$cell_sample = paste0(cell_names_all,'_',data_run_subset$sample )

  data_run_subset$sample_cell = paste0(data_run_subset$Sample,' ', cell_names_all)
  
  
  path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Cluster/PCA30/res3/Data/Labels.csv'
  data_run_subset = addOldLabels(path,data_run_subset, 'OldCellType')
  
  path = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/T Cell/Cluster/PCA30/res3.5//Data/labels.csv'
  data_run_subset = addOldLabels(path,data_run_subset, 'OldTcellLabels')

  path = "/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK_RemoveRiboFeatures/Cluster/PCA30/res3//Data/labels.csv"
  data_run_subset = addOldLabels(path,data_run_subset, 'OldNKLabels')
  
  path = "/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/Mono_DC/Cluster/PCA30/res1.6//Data/labels.csv"
  data_run_subset = addOldLabels(path,data_run_subset, 'OldMonoLabels')
  
  ##
 


  #data_run_subset = FindClusters(data_run_subset,resolution = resolution_val_subset)
  plot = DimPlot(data_run_subset,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
  print(plot)
  
  
  
  next
  path = paste0(filepath_cluster,'/Data/')
  dir.create(path,recursive = T)
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  all_cells = as.data.frame(colnames(data_run_subset_label))
  colnames(all_cells) = 'cell'
  all_cells$sample = data_run_subset_label$sample
  all_cells$sample_cell = data_run_subset_label$sample_cell
  all_cells$cell_type = Idents(data_run_subset_label)
  write.csv(all_cells, file = paste0(path,'labels','.csv'), quote=FALSE, row.names= F)
  
}

gene_name_vector= rownames(data_run_subset)
gene_name_vector[grep("HLA-DR", gene_name_vector)]
gene_list = 'LY6C1, LY6G6D, MMP9, TIMP1, LGALS3, S100A6, S100A11, ISG15, MS4A4C'
gene_list = c('HLA-A', 'HLA-B', 'HLA-C')
gene_list = unlist(strsplit(gene_list, ",")) 
gene_list = trimws(gene_list, which = c("both"), whitespace = "[ \t\r\n]")
gene_list [ !(gene_list %in% gene_name_vector )]


filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',PCA_dim_subset,'/res',resolution_val_subset,'/' )

data_run_subset = addMetaData(data_run_subset, metaData)
data_run_subset = load_Doublets(data_run_subset)
#data_run_subset = load_CellLabel(data_run_subset)
print(unique(data_run_subset$GeneralCellType))

data_run_subset$kit = data_run_subset$`10X kit`
data_run_subset$FeatureLessThan400 = data_run_subset$nFeature_RNA < 400

##################
data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
##################



##

plot = DimPlot(data_run_subset,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

plot = DimPlot(data_run_subset_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

if(celltype == 'T Cell'){
  TSCM_res4 = read.csv(paste0(filepath_cluster,'Data/','TSCM_res4','.csv'))

  new_Idents = as.character(Idents(data_run_subset_label))
  new_Idents[colnames(data_run_subset_label) %in% TSCM_res4$x] = 'TSCM'
  Idents(data_run_subset_label) = new_Idents
  
  remove_list = c('Remove','Pro Erythrocyte')
  data_run_subset_label = data_run_subset_label[,!(Idents(data_run_subset_label) %in% remove_list)]

}

unique(Idents(data_run_subset_label))
remove_list = c('Erythrocyte','DC/T-Cell DBL','Mono/T-Cell DBL','14','MK',
                'Remove','Pro Erythrocyte','Mono/CD8+ T Cell DBL','GMPC',
                'CD14+ Mono/T-cell DBL','CD14+ Mono/CD8+ T-cell DBL','DC/T-cell DBL','dDC','dMono','dMIP1a+ Mono',
                1:24,'dT-cell','CD16+ Mono/T-cell DBL','sDC','sMono')
remove_list = c('Remove','Erythrocyte')
#data_run_subset_label_remove = data_run_subset


data_run_subset_label_remove = data_run_subset_label[,!(Idents(data_run_subset_label) %in% remove_list)]

Idents(data_run_subset_label_remove) = as.character(Idents(data_run_subset_label_remove))
unique(Idents(data_run_subset_label_remove))
sort(table(Idents(data_run_subset_label_remove)))

plot = DimPlot(data_run_subset_label_remove,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)
# data_run_subset_label$Ident_32_33 = 'other'
# data_run_subset_label$Ident_32_33[Idents(data_run_subset_label) == '32'] = '32'
# data_run_subset_label$Ident_32_33[Idents(data_run_subset_label) == '33'] = '33'
# data_run_subset_label$Ident_32_33 = factor(data_run_subset_label$Ident_32_33, levels = c('other','32','33'))
# plot = DimPlot(data_run_subset_label,pt.size = 1, reduction = "umap",
#                label = F,group.by  = 'Ident_32_33')
# print(plot)

#
groupBy_list = c('sample','Diagnosis','kit','OldCellType',
                 'Treatment','Batch','LowCount',
                 'Doublet','FeatureLessThan400','Sample Type','Technology')
#groupBy_list = c('sample')
featurePlot_list = c('percent.mt','nCount_RNA','G2M.Score','S.Score')
splitBy_list = NA

filepath_cluster = paste0( folder_subcluster, 'Cluster/', 'PCA',PCA_dim_subset,'/res',resolution_val_subset,'/' )

#


plotAll(data_run_subset, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF =T, 
        clusterTF =F, markersTF = F, 
        groupBy = groupBy_list, splitBy = splitBy_list,featurePlot_list = featurePlot_list,
        PCA_dim = PCA_dim_subset,resolution_val = resolution_val_subset,
        pt.size = 2)

# data_run_subset_label = data_run_subset_label[,Idents(data_run_subset_label) != 14]


plotAll(data_run_subset_label, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = T, 
        clusterTF = F, markersTF = T,  
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim_subset,resolution_val = resolution_val_subset, 
        pt.size = 2,
        str = '_label')

plotAll(data_run_subset_label_remove, folder = folder_subcluster,
        sample_name,sampleParam = NA,
        cell_features = cell_features,
        label_TF = F,integrate_TF = F,  DE_perm_TF = T, 
        clusterTF = F, markersTF = T  ,
        groupBy = groupBy_list, splitBy = splitBy_list,
        PCA_dim = PCA_dim_subset,resolution_val = resolution_val_subset, 
        pt.size = 2,
        str = '_label_remove')

cell_features_file = paste0('/home/sujwary/Desktop/scRNA','/Data/Cell_IDS.xlsx')
cell_features = read_excel(cell_features_file)
cell_features = cell_features[cell_features$Plot_marker == 1,]

PlotKnownMarkers(data_run_subset, folder = paste0(filepath_cluster,'Cell Type/CellMarkers/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F, plotAll= F)

PlotKnownMarkers(data_run_subset, folder = paste0(filepath_cluster,'Cell Type/AllMarkers/'), 
                 cell_features = cell_features,
                 plotType ='FeaturePlotFix' , str = '',markerSize = 1, 
                 prefix_logFC = F, plotAll = T)


gene_list = c('AXL','GAS6')
gene_list = c('KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4', 'KIR2DL5','KIR2DL5A', 'KIR2DL5B',
              'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4DEL', 
              'KIR2DS4INS','KIR2DS5', 'KIR3DL1', 'KIR3DL2', 'KIR3DL3','KIR3DS1', 'KIR2DP1', 'KIR3DP1')


gene_list = c('TIGIT', 'LAG3', 'HAVCR2', 'B3GAT1', 'TNFSF10', 'PRF1', 'GZMK', 'GZMB', 'GZMH', 'ZBTB16', 'KLRC2')
gene_list = c('HLA-A', 'HLA-B', 'HLA-C','RELA')

gene_list = 'NKG7,PRF1, GZMB,GNLY, CCL5, FCGR3A, FCER1G,TYROBP,CX3CR1, ADGRG1,KLRD1, KLRF1,ITGAL, ZEB2,SRGN, GZMH, IFI16, RUNX3,IFITM2,IFITM1,GZMK,SELL,XCL1,XCL2,CD44,KLRC1,CD2,NCAM1,NR4A2,FOS,IER2,JUN,FOSL2,MAP3K8,NFKBIA,DUSP1,JUND,IRF1,CD3D,CD3G,GZMH,CD52,TRAC,TRGC2,LMNA,KDM6B,AREG,KDM6B,AREG,TNFAIP3,SLA2,OTULIN,BIRC2,TGFB1,REL,MAP3K8,NR4A2,MAFF,NFKBIA,CXCR4,SAT1,CD3D,CXCR4,AREG,CCL5,AREG,GZMK,XCL1,FOS,CCL3,NFKBIA,XCL2,CCL5, JUNB, FOSB,IER2,CD69,JUN,KLRB1,DUSP1,KLRF1,NR4A2,CCL4,JUNB,CD69,DUSP1,JUND,IER2,NR4A2,DUSP2,GZMK,JUN,JUNB,IER2,NKG7,DUSP1,FOSB,CCL3,CCL5,CFL1,FOS,CD69,KLRG1,KLRD1,KIR3DX1,KIRDL1'
FeaturePlot_GeneList(data_run_subset,gene_list,
                     folder = paste0(filepath_cluster,'Cell Type/12-10-20/'),
                     FeaturePlotFix = T,str = '', label = F)

remove_list = c('Remove','14','Erythrocyte','DC/T-Cell DBL',
                'Mono/CD8+ T Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL')
data_run_subset_label = data_run_subset_label[,!(Idents(data_run_subset_label) %in% remove_list)]
StackedVlnPlotHelper(data_run_subset_label,gene_list,
                     folder_heatMap = paste0(filepath_cluster,'Plots/Violin/'),
                     filename = paste0('AXL_GAS6_AllCluster.png'))


var = 'OldMonoLabels'
var = 'OldTcellLabels'
dir.create(paste0(filepath_cluster,'/Data/'))
cell_table = table( (data_run_subset@meta.data[[var]]),as.character(Idents(data_run_subset) ))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/',var,'_Table.csv'))
cell_table = table( data_run_subset$Technology,Idents(data_run_subset))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','TechnologyByIdent_Table.csv'))


cell_table = table( data_run_subset_label@meta.data[[var]],Idents(data_run_subset_label))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','OldLabels_Table_datalabel.csv'))
cell_table = table( data_run_subset_label$sample,Idents(data_run_subset_label))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','sampleByIdent_Table_datalabel.csv'))

cell_table = table( data_run_subset_label$Technology,Idents(data_run_subset_label))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','TechnologyByIdent_label.csv'))

cell_table = table( paste0(data_run_subset_label_remove$Diagnosis,' ',data_run_subset_label_remove$`Sample Type`,' ',data_run_subset_label_remove$Treatment),
                    Idents(data_run_subset_label_remove))
cell_table = cell_table/(rowSums(cell_table))
write.csv(cell_table, file = paste0(filepath_cluster,'/Data/','DiagnosisByIdent_prop_label.csv'))


####################
### Plot DE
####################
str = '_label'
str = ''
filepath = paste0(filepath_cluster
                  ,'Features',str,
                  '.csv')
DE = read.csv(file = filepath)
data_input = data_run_subset


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
  DE_cluster = DE_cluster[1:30,]
  
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
##

## Random plotting
group = 'OldCellType'
group = 'OldTcellLabels'
group = 'OldNKLabels'
group = 'OldMonoLabels'
#data_run_subset[[group]] = factor(as.character(data_run_subset[[group]]))
data_run_subset_remove = data_run_subset
data_run_subset_remove = data_run_subset[,data_run_subset[[group]] !='']


pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim_subset,'_res',resolution_val,
                                            '_GroupBy',group,'_label','.png'))
png(file=pathName,width=2000, height=1000)

plot = DimPlot(data_run_subset_remove, pt.size = 2, reduction = "umap",label = T,
               group.by  = group, label.size = 12)
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


### Plot DE

pathways.hallmark <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/h.all.v7.2.symbols.gmt')
pathways.c2 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c2.all.v7.2.symbols.gmt')
pathways.c5 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c5.all.v7.2.symbols.gmt')
pathways.c7 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c7.all.v7.2.symbols.gmt')
pathways.c8 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c8.all.v7.2.symbols.gmt')
pathways = c(pathways.hallmark,pathways.c2,pathways.c5, pathways.c7,pathways.c8)

str = ''

filepath = paste0(filepath_cluster
                  ,'Features',str,
                  '.csv')
DE = read.csv(file = filepath)
data_input = data_run_subset

DE = DE[DE$p_val_adj < 0.05,]
DE = DE[!grepl("MT-", DE$gene),]
DE = DE[!grepl("^RP[SL]", DE$gene),]

cluster_list = sort(unique(DE$cluster))

#cluster_list = cluster_list[4:length(cluster_list)]
for (cluster in cluster_list){
  DE_cluster = DE[DE$cluster == cluster,]
  ranks <- DE_cluster$avg_logFC
  names(ranks) <- DE_cluster$gene
  
  DE_cluster = DE_cluster[DE_cluster$avg_logFC > 0 ,]
  DE_cluster = DE_cluster[order(-DE_cluster$avg_logFC),]
  DE_cluster = DE_cluster[1:40,]
  
  DE_cluster = DE_cluster[rowSums(is.na(DE_cluster)) != ncol(DE_cluster),]
  gene_list = DE_cluster
  #gene_list = as.data.frame(as.character(gene_list[!is.na(gene_list)]))
  #colnames(gene_list) = 'Markers'
  
  
  gene_list$Cell = cluster
  gene_list$Plot_marker = 1
  gene_list$Markers = gene_list$gene
  path = paste0(filepath_cluster,'DE/Cluster_',cluster,'VsAll',str,'/')
  PlotKnownMarkers(data_input, 
                   folder = path, 
                   cell_features = gene_list,
                   plotType ='FeaturePlotFix' , 
                   prefix_logFC = T, str = '',markerSize = 1)
  next
  if (nrow(gene_list) > 3){
    path = paste0(filepath_cluster,'DE/GSEA/')
    dir.create(path,recursive = T)
    path = paste0(path,cluster,'VsAll',str,'GSEA.csv')
    
    
    
    next
    fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
    fgseaRes  = fgseaRes[order(fgseaRes$padj, fgseaRes$pval),]
    fgseaRes$leadingEdge =  sapply( fgseaRes$leadingEdge , paste0, collapse=",")
    write.csv(fgseaRes, file = path,row.names=F)
    
  }
  
  
}


### Plot DE permutations
data_input = data_run_subset_label
data_input = data_run_subset
cluster_list1 = sort(unique(Idents(data_input)))
cluster_list2 = cluster_list1

cluster_list1 = c('0','2','15')
cluster_list2 = c('12','13')

cluster1 = 11
cluster2 = 5
str = '_label_remove'
str = '_label'
str = ''
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

# Heatmap

newIdents = as.character(Idents(data_run_subset))

newIdents[newIdents == '8'] = '8_11'
newIdents[newIdents == '11'] = '8_11'
Idents(data_run_subset) = newIdents


gene_name_vector= rownames(data_run_subset)

gene_list = 'PTPRC, CD3D, CD3G, CD3E, CD4, CD8A, CD8B, IL2RA, FOXP3, CTLA4, 
  TIGIT, LAG3, SELL, CCR7, CD27, CD28, IL7R, NOSIP, LTB, CDK6, 
  LEF1, TCF7, BCL2, MYC, FCGR3A, CX3CR1, CXCR4, ADGRG1, GZMA, 
  GZMK, GZMH, GZMB, CCL3, CCL4, CCL5, PRF1, NKG7, RORA, KLRB1, 
  CD74, S100A4, S100A11, LGALS1, VIM, ANXA2, CD69, FAS, HNRNPLL, 
  GNLY, ITGB1, JUN, JUNB, JUND, FOS, FOSB, DUSP1, TNFAIP3, NFKBIA, 
  TSC22D3, CTSW, LDHB, ACTB, ACTG1, LIMS1, PASK, SOCS3, GATA3, 
  CCL4L2, CMC1, HLA-DRA, HLA-DRB1, HLA-DPB1, HLA-DPA1'
gene_list = unlist(strsplit(gene_list, ",")) 
gene_list = trimws(gene_list, which = c("both"), whitespace = "[ \t\r\n]")
gene_list [ !(gene_list %in% gene_name_vector )]

data_run_subset_label = ScaleData(data_run_subset_label,features =  rownames(data_run_subset_label))

unique(Idents(data_run_subset_label))

remove_list = c('MK','14','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL','Erythrocyte','DC/T-Cell DBL',
                '12','0','11','17','18','21','T/NK Doublet')
data_run_subset_label = data_run_subset_label[,!(Idents(data_run_subset_label) %in% remove_list)]

gene_list = gene_list = c('AXL','GAS6')
folder = paste0(filepath_cluster,'/Plots/Heatmap/')
dir.create(folder,recursive = T)
pathName <- paste0(folder,paste0('HeatMap', '_AXL_GAS6','.pdf'))
pdf(file=pathName, height=12, width=20)
#png(file=pathName)
#png(file=pathName,width=2000, height=1000,res = 100)

plot = DoHeatmap(data_run_subset_label, features = gene_list)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()


########
## DE
########

newIdents = as.character(Idents(data_run_subset))
newIdents[newIdents == '8'] = '8_11'
newIdents[newIdents == '11'] = '8_11'
Idents(data_run_subset) = newIdents

data_run_input = data_run_subset_label
# Seurat
ident1 = '8_11'
ident2 = '0'
ident1 = '1'
ident2 = 'cCD56dim'

Features = FindMarkers(data_run_input, ident.1 = ident1, ident.2 = ident2
                       ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE)

path = paste0(filepath_cluster,'DE/')
dir.create( path, recursive = TRUE)
path = paste0(path, 'DE ',ident1,' Vs ', ident2,'.csv')
print(path)
write.csv(Features, file = path,row.names=TRUE)

Features = read.csv(path)
str = ''
Features$Cell = ''
Features$Plot_marker = 1
Features$Markers = Features$X
Features = Features[order(Features$logFC),]

path = paste0(filepath_cluster,'DE/Cluster_',ident1,'Vs',ident2,str,'/')
PlotKnownMarkers(data_run_input, 
                 folder = path, 
                 cell_features = Features,
                 plotType ='FeaturePlotFix' , 
                 prefix_logFC = T, str = '',markerSize = 1)

plot = DoHeatmap(data_run_subset_label, features = gene_list)
plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24),
  text = element_text(size = 20))
print(plot)
dev.off()


# #######################
doubletMethod = 'Doublet3Methods'
mean_doublets = 100*sum((data_run_subset@meta.data[,doubletMethod]))/ncol(data_run_subset)
title_str = paste0('% Doublets: ',specify_decimal(mean_doublets,2))
pathName <- paste0(filepath_cluster,paste0(doubletMethod,'_GeneralCellType','.png'))
png(file=pathName,width=1000, height=500,res = 100)
plot = DimPlot(data_run_subset,
               pt.size = ifelse(data_run_subset@meta.data[,doubletMethod] == T, 1, 0.5), cols = c('navajowhite2','tomato2'),
               reduction = "umap",label = F, group.by  = doubletMethod) +
  ggtitle(title_str)
plot = plot +  DimPlot(data_run_subset,pt.size = 0,
                       reduction = "umap",label = F, group.by  = 'GeneralCellType')
print(plot)
dev.off()

## Barplot for CD3+/CD14+ cluster by sample

cells  = colnames(data_run_subset)
cells = cells[Idents(data_run_subset) %in% c(8,17)]
samples = data_run_subset$sample[Idents(data_run_subset) %in% c(8,17)]

TcellMono_cells = data.frame(cbind(cells,samples))
TcellMono_cells$cell_sample = data_run_subset$cell_sample[Idents(data_run_subset) %in% c(8,17)]
write.csv(TcellMono_cells, file = paste0(filepath_cluster,'/Stats/','TcellMono_cells'))


# Matrix for NMF
cell_list = c("Naive CD8+ T-cell","cTreg" ,"GZMH+ GZMB+ CD8+ T-cell","CD8+ TCM",
              "GZMK+ CCL3+ CCL4+ CD8+ T-cell", "IFN+ CD4+ T-cell" ,
"TSCM"    ,"Naive CD4+ T-cell"  ,"eTreg"  , "CD4+ TCM" ,"Th2",  "TEMRA" ,                   
"GZMK+ CD8+ T-cell" ,"Stim Naive CD4+ T-cell" , "CCL5+ CD4+ T-cell"  ,"aTh17" , "Th17","TRM")                
celltype = 'TCell_subset'

cell_list = c("cDC2" ,                     "sDC"   ,                  
"prDC"     ,                "cDC1"  , 
"CD16+ Mono" ,               "TGFb1+ CD14+ Mono",      
"SELL+ CD14+ Mono" ,          "CD14+ CD16+ Mono"   ,     
"IFN+ Mono"      ,           "MIP1a+ CD14+ Mono" ,'TIMP+ CD14+ Mono' )
celltype = 'MonoDC_subset'

cell_list = c("Tgd"  ,"aCX3CR1+ GZMB+ CD56dim" ,   "aCXCR4+ CD56dim"    
              ,"aGZMK+ CCL3+ CD56dim",  "cCD56dim" ,"CD56bright"   ,"NFkB-high",'aCCL3+ GZMK+')  
celltype = 'NK_subset'


unique(Idents(data_run_subset_label))


cell_list[!(cell_list %in% unique(Idents(data_run_subset_label)))]

unique(Idents(data_run_subset_label))[!(unique(Idents(data_run_subset_label)) %in%cell_list )]

data_input = data_run_subset_label[,Idents(data_run_subset_label) %in% cell_list]
plot = DimPlot(data_input,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

data_matrix = data_input@assays[["RNA"]]@data
data_matrix_var = data_matrix[rownames(data_matrix) %in% data_input@assays[["RNA"]]@var.features,]
path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_subset','.tsv')
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')
path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_varFeature_subset','.tsv')
write.table(data_matrix_var, 
            file=path, 
            quote=FALSE, sep='\t')

data_input = data_run_subset
data_matrix = data_input@assays[["RNA"]]@data
data_matrix_var = data_matrix[rownames(data_matrix) %in% data_input@assays[["RNA"]]@var.features,]
path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'.tsv')
write.table(data_matrix, 
            file=path, 
            quote=FALSE, sep='\t')
path = paste0(filepath_cluster,'Data/matrix_',
              celltype,'_varFeature','.tsv')
write.table(data_matrix_var, 
            file=path, 
            quote=FALSE, sep='\t')

# UmapCoord and metadata

sample = as.character(data_run_subset$sample)
ident = as.character(Idents(data_run_subset))
UmapCoord = data_run_subset@reductions[["umap"]]@cell.embeddings
cell =  colnames(data_run_subset)
label = as.character(Idents(data_run_subset_label))
old_ident = as.character(data_run_subset$OldCellType)
old_ident_T = as.character(data_run_subset$OldTcellLabels)
old_ident_NK = as.character(data_run_subset$OldNKLabels)
old_ident_Mono = as.character(data_run_subset$OldMonoLabels)


output = cbind(cell,sample,ident,label,UmapCoord,old_ident, old_ident_T,old_ident_NK,old_ident_Mono)
#output = output[,c('sample')]
#output = data_run_subset_label@meta.data
dir.create(paste0(filepath_cluster,'Data/'))
path = paste0(filepath_cluster,'Data/data.csv')
write.csv(output, file = path,row.names=TRUE,col.names=TRUE)



## AUCell

exprMatrix = data_run_subset@assays[["RNA"]]@counts

library(AUCell)
library(GSEABase)
gmtFile <- paste0('/home/sujwary/Desktop/scRNA/Data/msigdb.v7.1.symbols.gmt')
geneSets <- getGmt(gmtFile)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE) 

cells_assignment

path = paste0(folder_subcluster,'cells_AUC','.Robj')
save(cells_AUC,file= path)

cellAssigned <- cells_assignment$Oligodendrocyte_Cahoy$assignment

geneSetName <- rownames(cells_AUC)[grep("Oligodendrocyte_Cahoy", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)

assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

miniAssigMat <- assignmentMat[,sample(1:ncol(assignmentMat),100)]
library(NMF)
aheatmap(miniAssigMat, scale="none", color="black", legend=FALSE)


## scVelo

data_input = data_run_subset
#data_input = data_input[,!(Idents(data_input) %in% c(0,8,11,17,18,21,22))]
data_input$sample_cell = paste0(data_input$Sample,'_', colnames(data_input))

dir.create(paste0(filepath_cluster,'Data/'))

write.csv(unname(data_input$sample_cell), 
          file = paste0(filepath_cluster,"Data/cellID_obs.csv"))
write.csv(Embeddings(data_input, reduction = "umap"), 
          file = paste0(filepath_cluster,"Data/cell_embeddings.csv"))
write.csv(Idents(data_input), file = paste0(filepath_cluster,"Data/clusters.csv"))

library(loomR)

for(j in 1:ncol(data_input@meta.data)){
  if(is.factor(data_input@meta.data[,j]) == T){
    data_input@meta.data[,j] = as.character(data_input@meta.data[,j]) # Force the variable type to be character
    data_input@meta.data[,j][is.na(data_input@meta.data[,j])] <- "NA"
  }
  if(is.character(data_input@meta.data[,j]) == T){
    data_input@meta.data[,j][is.na(data_input@meta.data[,j])] <- "N.A"
  }
}
meta_data = data_input@meta.data
data_input@meta.data = meta_data[,c("orig.ident","nCount_RNA" ,"percent.mt" , "seurat_clusters")]
pfile <- as.loom(data_input, filename = paste0(filepath_cluster,"Data/data.loom"), verbose = FALSE)
 
# Proportion analysis

