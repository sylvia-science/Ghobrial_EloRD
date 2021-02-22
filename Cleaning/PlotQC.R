library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
output_folder = '/disk2/Projects/EloRD_Nivo_Oksana/Output/Soup_MT_C100/'
output_folder = '/disk2/Projects/MMRF/Output/QC/'

data_folder = "/home/sujwary/Desktop/scRNA/Data/"
filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'

data_folder = '/disk2/Projects/EloRD_Nivo_PBMC/Data/'
filename_metaData = '/disk2/Projects/EloRD_Nivo_PBMC/MetaData/metaData_EloRD_Nivo_PBMC.xlsx'

data_folder = '/disk2/Projects/EloRD_Nivo_Oksana/Data/'
filename_metaData = '/disk2/Projects/EloRD_Nivo_Oksana/MetaData/10X Sequenced Samples.xlsx'

data_folder = '/disk2/Projects/MMRF/Data/'
filename_metaData = '/disk2/Projects/MMRF/MetaData/MetaData.csv'


metaData = read_excel(filename_metaData)
metaData = read.csv(filename_metaData)

metaData = metaData[metaData$Run== 1,]
#metaData = metaData[metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]
#metaData = metaData[metaData$Study == 'Nivo',]
#metaData = metaData[metaData$`10X kit` == 'Microwell-seq',]



#filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample','_parameters.xlsx')
#sampleParam <- read_excel(filename_sampleParam)

#sampleParam = sampleParam[sampleParam$Sample %in% metaData$Sample,]

i = 1

# Soup + MT + Normal threshold

sample_list = c('GL1497BM', 'GL1160BM', 'GL2923BM', 'GL3404BM', 'NBM6CD138N', 'NBM12CD138N', 'GL2185BM', 'GL3417BM', 'GL2653BM')

#sample_list = c('GL3404BM')
i = 11
run = T
for (i in 1:nrow(metaData)){
  sample_name = metaData$Sample[i]

  path = paste0(folder,'/',sample_name,'.Robj')
  if (file.exists(path)){
    print('File Exists')
    #next
  }
  
  #sample_name = metaData$Sample[i]
  #sample_name = sample_list[i]
  #sample_name = 'GL1160BM'
  print(sample_name)
  #percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  #RNA_features_min = sampleParam$RNA_features_min[sampleParam['Sample'] == sample_name]
  #RNA_features_max = sampleParam$RNA_features_max[sampleParam['Sample'] == sample_name]
  
  #filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  file_list = list.files(path = data_folder)
  HCL_list = c('Adult-Bone-Marrow1','Adult-Bone-Marrow2',
               'Adult-Peripheral-Blood1','Adult-Peripheral-Blood2','Adult-Peripheral-Blood3','Adult-Peripheral-Blood4')
  if (sample_name %in% HCL_list){
    filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_dge.txt",sep = "")
    
    set  = ' '
    data_i_raw = read.table(file = filename,row.names = 1,header = T, sep = '')
    if (ncol(data_i_raw) == 1){
      data_i_raw = read.table(file = filename,row.names = 1,header = T, sep = ',')
    }
    nrow(data_i_raw)
    ncol(data_i_raw)
    data_i_raw = CreateSeuratObject(data_i_raw,  project = "BM",min.cells = 3, min.features = 1)
    
  }else{
    filename = file_list[startsWith(file_list,paste0(sample_name,'_raw') )]
    path = paste0(data_folder,filename)
    
    file.exists(path)
    data_i_raw = Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
    data_i_raw = CreateSeuratObject(counts = data_i_raw, project = "BM", min.cells = 3, min.features = 1)
  }
  colSum_list = colSums(data_i_raw ) # Needs to be from Matrix library
  mincount = 100
  keep = colSum_list >= 100
  print(sum(keep))
  print(sum(!keep))
  if (sum(!keep) < 3){
    mincount = 200
    
  }
  
  
  data_i_raw = data_i_raw[!is.na(rownames(data_i_raw)),]
  
  data_i_raw[["percent.mt"]] <- PercentageFeatureSet(data_i_raw, pattern = "^MT-")
  data_i_raw = data_i_raw[,!is.na(data_i_raw$percent.mt)]
  max(data_i_raw$percent.mt)
  
  folder = paste0(output_folder,'/')
  dir.create(folder,recursive = T)
  
  pathName <- paste0(folder, sample_name,'_QC_PreFilter.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data_i_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0)
  print(plot)
  dev.off()
  
  countSum_min = min(colSum_list)
  #next
  colSum_list = colSums(data_i_raw )
  keep = colSum_list >=mincount
  data_i_filtered = data_i_raw[,keep]
  
  #data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  


  
  folder = paste0(output_folder,'/')
  pathName <- paste0(folder, sample_name, '_QC_PostFilter.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data_i_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0,split.by = NULL)
  print(plot)
  dev.off()
  
  data_i_filtered = data_i_filtered[, data_i_filtered$percent.mt < 15]
  
  folder = paste0(output_folder,'/')
  pathName <- paste0(folder,sample_name, '_QC_PostFilterMT.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data_i_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0,split.by = NULL)
  print(plot)
  dev.off()
  
  
  
}
