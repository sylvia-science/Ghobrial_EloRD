library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')
# Run EmptyDrops on individual samples

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
filename_metaData = '/disk2/Projects/EloRD_Nivo_PBMC/MetaData/metaData_EloRD_Nivo_PBMC.xlsx'

data_folder = '/disk2/Projects/EloRD_Nivo_Oksana/Data/'
output_folder = '/disk2/Projects/EloRD_Nivo_Oksana/Output/EmptyCells/'
filename_metaData = '/disk2/Projects/EloRD_Nivo_Oksana/MetaData/10X Sequenced Samples.xlsx'

data_folder = '/disk2/Projects/EloRD_Nivo_Oksana/Data/'
output_folder = '/disk2/Data/EmptyCells/'
filename_metaData = '/disk2/Projects/MMRF/MetaData/MetaData.csv'


metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

#metaData = metaData[metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[rowSums(is.na(metaData)) != ncol(metaData), ]
#metaData = metaData[metaData$`10X kit` == 'Microwell-seq',]
#metaData = metaData[metaData$Study == 'Nivo',]

sample_name = 'GL2923BM'
#filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_parameters.xlsx')
#sampleParam <- read_excel(filename_sampleParam)
i = 29

for (i in 1:nrow(metaData) ){

  sample_name = metaData$Sample[i]
  print(sample_name)
  
  file = paste0(output_folder,sample_name,'/', 'emptyDrops','.csv')
  if(file.exists( file)){
    print('File already exists')
    next
  }

  
  HCL_list = c('Adult-Bone-Marrow1','Adult-Bone-Marrow2',
               'Adult-Peripheral-Blood1','Adult-Peripheral-Blood2','Adult-Peripheral-Blood3','Adult-Peripheral-Blood4')
  # Files from the HCL need to be read in differently
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
    file_list  = list.files(path = data_folder)
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
  
  countSum_min = min(colSum_list)
  #next
  keep = colSum_list >=mincount
  data_i_filtered = data_i_raw[,keep]
  
  # If data was already filtered,just says all are real cells
  if (countSum_min > 400){
    rownames = colnames(data_i_filtered)
    br_e = data.frame(rownames)
    br_e$is_cell <- T
    br_e$ncount = data_i_filtered$nCount_RNA
    br_e$LogProb = 0
    br_e$FDR  = 1
    
    folder = paste0(output_folder,sample_name,'/')
    dir.create(folder, recursive =T)
    write.csv(br_e, file = paste0(folder,'emptyDrops','.csv'),row.names = FALSE)
    
    print('Min count sum is larger than 400')
    next
  }

  
  ## Empty drops
  counts= data_i_raw@assays[["RNA"]]@counts
  br.out = barcodeRanks(counts,lower=mincount)
  e.out = emptyDrops(counts,lower=mincount) # Outputs NAs for cell with lower than 100
  
  
  br = data.frame(br.out@rownames,br.out$rank, br.out$total)
  e = data.frame(e.out@rownames, e.out$FDR,e.out$LogProb)
  colnames(br) = c('rownames','rank','total')
  colnames(e) = c('rownames','FDR','LogProb')
  
  br_e <- merge(br,e,by="rownames")
  br_e$is_cell <- br_e$FDR <= 0.01
  br_e$ncount = data_i_raw$nCount_RNA

  #ggplot(br_e) + geom_point(aes(rank,total,color=is_cell)) + scale_y_log10() + scale_x_log10()
  
  
  br_e_sorted = br_e[match(as.character(rownames(data_i_filtered@meta.data)), as.character(br_e$rownames)   ),]
  
  folder = paste0(output_folder,sample_name,'/')
  dir.create(folder, recursive =T)
  write.csv(br_e, file = paste0(folder,'emptyDrops','.csv'),row.names = FALSE)
  
  next
  #######################################################
  ## Run pipeline before and after removing empty cells
  #######################################################
  data_i_filtered$emptyProb = br_e_sorted$LogProb
  data_i_filtered$is_cell = br_e_sorted$is_cell
  data_i_filtered[["percent.mt"]] <- PercentageFeatureSet(data_i_filtered, pattern = "^MT-")
  
  
  data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
  data_matrix_raw = data_i_raw@assays[["RNA"]]@counts
  cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))
  
  pt.size = 0.8
  
  print('Plot')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_is_cell','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'is_cell'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_emptyProb','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'emptyProb'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_run,pt.size = pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0(output_folder,sample_name,'/','Featureplots/PreEmpty/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_run, feature = gene, folder =folder,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    file_str = ''
    pathName = paste0(folder,gene,file_str,'.png')
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  
  data_i_filtered_isCell = data_i_filtered
  data_i_filtered_isCell = data_i_filtered_isCell[,data_i_filtered_isCell$is_cell  == T]
  
  data_i_filtered_isCell = NormalizeData(data_i_filtered_isCell, normalization.method = "LogNormalize", scale.factor = 10000)
  data_i_filtered_isCell = FindVariableFeatures(data_i_filtered_isCell, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_isCell = ScaleData(data_i_filtered_isCell)
  data_i_filtered_isCell = RunPCA(data_i_filtered_isCell,npcs = 30)
  data_i_filtered_isCell = FindNeighbors(data_i_filtered_isCell, dims = 1:30)
  data_i_filtered_isCell = FindClusters(data_i_filtered_isCell)
  data_i_filtered_isCell = RunUMAP(data_i_filtered_isCell, dims = 1:30)
  
  print('Plot')
  folder = paste0(output_folder,sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PostEmpty_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_isCell,pt.size = pt.size, reduction = "umap",label = FALSE))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PostEmpty_Umap','_percent.mt','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_isCell,pt.size = pt.size, reduction = "umap", features = 'percent.mt'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PostEmpty_Umap','_is_cell','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_isCell,pt.size = pt.size, reduction = "umap", features = 'is_cell'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PostEmpty_Umap','_emptyProb','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_isCell,pt.size = pt.size, reduction = "umap", features = 'emptyProb'))
  dev.off()
  
  pathName = paste0(folder,sample_name,'_PostEmpty_Umap','_nFeature_RNA','.png')
  png(file=pathName,width=1000, height=1000)
  print(FeaturePlot(data_i_filtered_isCell,pt.size = pt.size, reduction = "umap", features = 'nFeature_RNA'))
  dev.off()
  
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0(output_folder,sample_name,'/','Featureplots/PostEmpty/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_isCell, feature = gene, folder =folder,
                          str = '',split = F, markerSize = 3,gene_TF = TRUE,title = '',saveTF = FALSE)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=24 ),
      axis.title.y = element_text(color="black", size=24),
      axis.text= element_text(color="black", size=24),
      legend.text=element_text(size=24),
      legend.title=element_text(size=24),
      text = element_text(size = 20)
    )
    
    file_str = ''
    pathName = paste0(folder,gene,file_str,'.png')
    pathName = paste0(folder,gene,'','.png')
    png(filename = pathName,width=2000, height=2000)
    print(plot)
    dev.off()
    remove(plot)
    
  }
  
  
}

# End main run
################################################
## Apply emptyDrops to Soup corrected


for (i in 1:nrow(metaData) ){

  sample_name = metaData$Sample[i]
  print(sample_name)
  Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  print(Scrublet_threshold)
  percent_mt = sampleParam$percent_mt_min[sampleParam['Sample'] == sample_name]
  
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/Soup/',sample_name,'/')
  path = paste0(folder,'/',sample_name,'_Soup_adjusted.Robj')
  data_soup = loadRData(path)
  
  
  row_val = rownames(data_soup)
  col_val = colnames(data_soup)
  
  ## Empty drops
  

  folder = paste0(output_folder,sample_name,'/')
  file_name = paste0(folder,'emptyDrops','.csv')
  br_e = read.csv(file_name)
  br_e_sorted = br_e[match(as.character(rownames(data_soup@meta.data)), as.character(br_e$rownames)   ),]
  
  
  data_soup$emptyProb = br_e_sorted$LogProb
  data_soup$is_cell = br_e_sorted$is_cell
  data_soup[["percent.mt"]] = PercentageFeatureSet(data_soup, pattern = "^MT-")
  
  data_soup_filtered = data_soup[,data_soup$percent.mt < 15]
  
  empty_preFilter =colnames(data_soup)[(data_soup$is_cell == F | is.na(data_soup$is_cell ))]
  empty_postFilter = colnames(data_soup_filtered)[(data_soup_filtered$is_cell == F | is.na(data_soup_filtered$is_cell ))]
  
  nCount_preFilter = data_soup$nCount_RNA[empty_preFilter]
  MT_preFilter = data_soup$percent.mt[empty_preFilter]
  MT_notempty_preFilter = data_soup$percent.mt[!(colnames(data_soup) %in% empty_preFilter)]
  
  
  nCount_postFilter = data_soup_filtered$nCount_RNA[empty_postFilter]
  MT_postFilter = data_soup_filtered$percent.mt[empty_postFilter]
  MT_notempty_postFilter = data_soup_filtered$percent.mt[!(colnames(data_soup_filtered) %in% empty_postFilter)]
  
  
  data_barplot_preFilter = data.frame(nCount_preFilter,MT_preFilter,'Pre MT')
  colnames(data_barplot_preFilter) = c('nCount','MT','Filter')

  if (length(nCount_postFilter) == 0){
    next
  }
  
  data_barplot_postFilter = data.frame(nCount_postFilter,MT_postFilter,'Post MT')
  colnames(data_barplot_postFilter) = c('nCount','MT','Filter')
  data_barplot = rbind(data_barplot_preFilter,data_barplot_postFilter)
  
  ##
  
  data_barplot_notEmpty_pre = data.frame(MT_notempty_preFilter,'Pre MT')
  colnames(data_barplot_notEmpty_pre) = c('MT','Filter')
  
  
  data_barplot_notEmpty_post = data.frame(MT_notempty_postFilter,'Post MT')
  colnames(data_barplot_notEmpty_post) = c('MT','Filter')
  
  data_barplot_notEmpty = rbind(data_barplot_notEmpty_pre,data_barplot_notEmpty_post)
  
  ##

  subtitle_str = paste0('Num Empty Cells Pre MT: ',length(empty_preFilter) ,'\n',
                        'Num Empty Cells Post MT: ',length(empty_postFilter))
  
  subtitle_str_notEmpty = paste0('Num Full Cells Pre MT: ',length(MT_notempty_preFilter) ,'\n',
                        'Num Full Cells Post MT: ',length(MT_notempty_postFilter))
  
  pathName = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/BarPlots/',sample_name,'_empty_counts','','.png')
  png(file=pathName,width=500, height=500)
  plot_count <- ggplot(data_barplot, aes(x=Filter, y=nCount)) + 
    geom_boxplot() + ggtitle('NCount in Empty Cells',  subtitle = subtitle_str) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #+
    #
  print(plot_count)
  dev.off()
  
  pathName = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/BarPlots/',sample_name,'_empty_MT','','.png')
  png(file=pathName,width=500, height=500)
  plot_MT <- ggplot(data_barplot, aes(x=Filter, y=MT)) + 
    geom_boxplot() + ggtitle('MT in Empty Cells',  subtitle = subtitle_str) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #+
  print(plot_MT)
  dev.off()
  
  pathName = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/BarPlots/',sample_name,'_full_MT','','.png')
  png(file=pathName,width=500, height=500)
  plot_notEmpty <- ggplot(data_barplot_notEmpty, aes(x=Filter, y=MT)) + 
    geom_boxplot() + ggtitle('MT in Full Cells',  subtitle = subtitle_str_notEmpty) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) #+
  print(plot_notEmpty)
  dev.off()

  
  
  

  
  
}
