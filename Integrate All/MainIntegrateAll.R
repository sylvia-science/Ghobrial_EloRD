# source('/home/sujwary/Desktop/scRNA/Code/Integrate All/MainIntegrateAll.R')
# Libraries
gc()

library(dplyr)
library(Seurat)
library(h5)
library(readxl)
library(ggplot2)
library(ggrepel)
library(stringr)
library(data.table)
library(hdf5r)
require(gridExtra)
#library(biomaRt)
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')


folder_base_input <- '/home/sujwary/Desktop/scRNA/Output/'
folder_base = '/home/sujwary/Desktop/scRNA/Output/'

filename_sampleParam <- '/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx'
sampleParam <- read_excel(filename_sampleParam)

filename_sampleParam_integrate <-'/home/sujwary/Desktop/scRNA/Data/sample_Combine_parameters.xlsx'
sampleParam_integrate = read_excel(filename_sampleParam_integrate)

filename_metaData <- '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData <- read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]
sample_type = 'PrePostEOTNBM_MT15'
integrate_merge = 'Merge'
folder_base_output = paste0('/home/sujwary/Desktop/scRNA/Output/',integrate_merge ,' All/',sample_type,'/')
dir.create( folder_base_output, recursive = TRUE)

print('Start')

filter = F
if (filter == FALSE){
  filterTF = '_filterF'
}else{
  filterTF = ''
}

data_list = vector(mode = "list", length = 0)
data_list = vector(length = 0)
sample_name_list = vector(mode = "list", length = 0)

for(i in 1:nrow(metaData)){
  
  sample_name = metaData$Sample[i]
  print(sample_name)
  Treatment = metaData$Treatment[i]
  dexa = metaData$'Dexa or not'[i]
  nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                            ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
  percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  if (metaData$`10X kit`[i] == 'Microwell-seq'){
    filename <- paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_dge.txt",sep = "")
    
    data_i = read.table(file = filename,row.names = 1,header = T)
    data_i = CreateSeuratObject(data_i)
  }else{
    filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
    filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_filtered_feature_bc_matrix.h5",sep = "")
    
    data_i = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
    
    data_i = CreateSeuratObject(counts = data_i, project = "BM", min.cells = 3,min.features = 0)
    
  }
  print(min(data_i$nCount_RNA))
  data_i[["percent.mt"]] = PercentageFeatureSet(data_i, pattern = "^MT-")
  
  
  if (filter){
    nFeature_RNA_min = as.numeric(nFeature_RNA_list[1])
    nFeature_RNA_max = as.numeric(nFeature_RNA_list[2])
    
    # Threshohold
    expr <- FetchData(object = data_i, vars = 'nFeature_RNA')
    data_i = data_i[, which(x = expr > nFeature_RNA_min & expr < nFeature_RNA_max)]
    expr <- FetchData(object = data_i, vars = 'percent.mt')
    data_i = data_i[, which(x = expr < percent_mt)]
  }
  
  if (integrate_merge == "Integrate"){
    data_i = NormalizeData(data_i, normalization.method = "LogNormalize", scale.factor = 10000)
    data_i = FindVariableFeatures(data_i, verbose = FALSE, selection.method = "vst", nfeatures =2000)
  }

  if (length(data_i$orig.ident) > 0){
    
    data_i@meta.data$orig.ident = paste0('data_',Treatment)
    data_i@meta.data$data_pre_renamed = data_i@active.ident
  
    data_i@meta.data$sample_name = sample_name
    data_i@meta.data$dexa = dexa
    data_i@meta.data$CHIP =  metaData$'CHIP or not'[i]
    
    data_i@meta.data$'Patient Number' = metaData$'Patient Number'[i]
    data_i@meta.data$Response = metaData$'Response'[i]
    data_i@meta.data$'10X kit' = metaData$'10X kit'[i]
    data_i@meta.data$'MM' = TRUE
    data_i@meta.data$'batch' = metaData$'batch'[i]
    data_i@meta.data$'SampleType' = metaData$'Sample Type'[i]
    
    print(data_i)
    
   
    if (metaData$`10X kit`[i] == 'Microwell-seq'){
      #browser()
      nFeature_RNA_list <- list(sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
                                ,sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name])
      percent_mt <- sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
      cluster_IDs = sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name] 
      
      resolution_val<- sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
      PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
      folder_base = '/home/sujwary/Desktop/scRNA/Output/'
      folder = makeFolders(folder_base,sample_name,filter = T,regress_TF = F, makeFolder_TF= F,nFeature_RNA_list,percent_mt)      
      print(paste0('folder: ', folder))
      
      data = loadRData(file=paste0(folder,paste0('data_run_PCAdim',PCA_dim,'.Robj')))
      data = FindClusters(data, resolution = resolution_val)
      data_run_label = label_cells(data,cluster_IDs)
      ind = which(colnames(data_i) %in% colnames(data_run_label))
      data_i@meta.data$'Cell_Label' = Idents(data_run_label) 
    }else{
      data_i@meta.data$'Cell_Label' = ''
    }
    
    #data_matrix = data_i@assays[["RNA"]]@counts
    #path_merge = paste0(folder_base_output,'PrePostEOTNBM_MT15_',sample_name,'.tsv')
    #write.table(data_matrix, file=path_merge, quote=FALSE, sep='\t')
    
    data_list = c(data_list, data_i)
    sample_name_list= c(sample_name_list,sample_name)
  }
      
}


for (i in 1:length(data_list_PCA)){
  data_i = data_list_PCA[i][[1]]
  print(sum(is.na(data_i@assays[["RNA"]]@data)))
}

write.csv(sample_name_list, file =paste0(folder_base_output,'sample_name_list','_',sample_type,'.csv'),row.names = FALSE,col.names = FALSE)


if (integrate_merge == 'Integrate'){
  
   for (i in 1:length(data_list_PCA)){
     data_i = data_list_PCA[i][[1]]
     print(length(data_i@assays[["RNA"]]@var.features))
   }
  
   
  
   # Regular Integration
   # dim:	Which dimensions to use from the CCA to specify the neighbor search space
   anchors = FindIntegrationAnchors(object.list = data_list, 
                                    dims = 1:30, k.filter =30, k.score = 30,
                                    anchor.features = 2000)
   path = paste0(folder_base_output,'anchors_',integrate_merge,'_',sample_type,filterTF,'.Robj')
   save(anchors,file=path)
   anchors = loadRData(path)
   
   # dim: Number of PCs to use in the weighting procedure
   data_integrate = IntegrateData(anchorset = anchors, 
                                  dims = 1:30,k.weight  = 100)
   path_integrate = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,filterTF,'.Robj')
   save(data_integrate,file=path_integrate)
   
   # RPCA
   features <- SelectIntegrationFeatures(object.list = data_list)
   data_list_PCA <- lapply(X = data_list, FUN = function(x) {
     print(unique(x$sample_name))
     x <- ScaleData(x, features = features, verbose = FALSE)
     x <- RunPCA(x, features = features,npcs = 30, verbose = T)
   })
   
  anchors_rpca = FindIntegrationAnchors(object.list = data_list_PCA,
                                        dims = 1:30, k.filter =30, k.score = 30,
                                        anchor.features = 2000,
                                        reduction = "rpca",scale = T)
  save(anchors_rpca,file=paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj'))
  anchors_rpca = loadRData(paste0(folder_base_output,'anchors_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj'))
  data_integrate_rpca = IntegrateData(anchorset = anchors_rpca,
                                      dims = 1:30,k.weight  = 100)

  
  path = paste0(folder_base_output,'data_',integrate_merge,'_rpca_',sample_type,filterTF,'.Robj')
  save(data_integrate_rpca,file=path)
  browser()
  # online fix
  # data_list <- lapply(X = data_list, FUN = function(x) {
  #   x <- SCTransform(x)
  # })
  # 
  # features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)
  # data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)
  # data_list <- lapply(X = data_list, FUN = RunPCA, verbose = FALSE, features = features)
  # anchors <- FindIntegrationAnchors(object.list = hca.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
  # data_integrate <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  # path  = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,'.Robj')
  # save(data_integrate,file=path)
  # 
}

browser()
if (integrate_merge == 'Merge'){
  data_merge = merge(x =  data_list[[1]],y = data_list[[2]], merge.data = FALSE)
  for (i in 3:length(data_list)){
    print(i)
    data_merge = merge(x =  data_merge,y = data_list[[i]],merge.data = FALSE)
  }
  
  path_merge = paste0(folder_base_output,'data_',integrate_merge,'_',sample_type,filterTF,'filtered','.Robj')
  
  save(data_merge,file=path_merge)
}

data_matrix = data_merge@assays[["RNA"]]@counts
path_merge = paste0(folder_base_output,'PrePostEOTNBM_MT15_All_AllFeatures_.tsv')
write.table(data_matrix, file=path_merge, quote=FALSE, sep='\t')


