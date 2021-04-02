library(Matrix)
library(readxl)
library(Seurat)
library(readxl)
library(DropletUtils)
library(ggplot2)
library(SoupX)

load_emptyDrops <- function(data,base) {
  data$emptyProb = NA
  data$is_cell = NA
  sample_list = unique(data$sample)
  for (i in 1:length(sample_list) ){
    
    sample_name = sample_list[i]
    print(sample_name)
    #browser()

    cell_names_all = colnames(data[,data$sample == sample_name])
    cell_names_all = sub("_.*", "", cell_names_all)
    

    
    folder = paste0(base,sample_name,'/')
    file_name = paste0(folder,'emptyDrops','.csv')
    if(!file.exists(file_name)){
      print(file_name)
      print('Does not exist')
      next
    }
    br_e = read.csv(file_name)
    br_e = br_e[!is.na(br_e$is_cell),]

    cell_names_sample = as.character(br_e$rownames)
    cell_names_sample = sub("_.*", "", cell_names_sample)
    br_e = br_e[cell_names_sample %in% cell_names_all,]
    cell_names_sample = cell_names_sample[cell_names_sample %in% cell_names_all]
    
    
    br_e_sorted = br_e[match(cell_names_sample,cell_names_all   ),]
    
    if (any(is.na(br_e_sorted$is.cell))){
      browser()
    }
    if (nrow(br_e_sorted) == 0){
      print('All NA values for empty cells')
      data$emptyProb[data$sample == sample_name] = 0
      data$is_cell[data$sample == sample_name] = T
      next
    }
    data$emptyProb[data$sample == sample_name] = br_e_sorted$LogProb
    data$is_cell[data$sample == sample_name] = br_e_sorted$is_cell

  }
  return (data)
}


load_Doublets <- function(data, folder) {
  data$scrublet = F
  data$doublet_finder = F
  data$scran_doublet = F
  data$scds_doublet = F
  data$Doublet = F
  data$Doublet2Methods = F
  data$Doublet3Methods = F
  data$doublet_scores = NA
  data$predicted_doublet_scrub = F
  sample_list = unique(data$sample)
  for (i in 1:length(sample_list) ){
    #browser()
    sample_name = sample_list[i]
    print(sample_name)
    #browser()
    #threshold=  sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
    path = paste0(folder,sample_name,'/doublet_summary','.csv')
    #browser()
    if (!file.exists(path)){
      print('No Doublet file found')
      next
    }
    doublet = read.csv(path,header = T)
    #browser()

    cell_names_all = colnames(data[,data$sample == sample_name])
    cell_names_all = sub("_.*", "", cell_names_all)
    
    cell_names_sample = as.character(doublet$Cell)
    cell_names_sample = sub("_.*", "", cell_names_sample)
    doublet = doublet[cell_names_sample %in% cell_names_all,]
    cell_names_sample = cell_names_sample[cell_names_sample %in% cell_names_all]
    
    doublet = doublet[match(cell_names_sample,cell_names_all),]
    
    #browser()
    #browser()
    data$scrublet[data$sample == sample_name] = doublet$scrublet
    data$doublet_finder[data$sample == sample_name] = doublet$doublet_finder
    data$scran_doublet[data$sample == sample_name] = doublet$scran_doublet
    data$scds_doublet[data$sample == sample_name] = doublet$scds_doublet
    data$Doublet[data$sample == sample_name]  = as.logical(as.character(toupper(doublet$Summary >= 2)))
    data$Doublet2Methods[data$sample == sample_name]  = as.logical(as.character(toupper(doublet$Summary >= 2)))
    data$Doublet3Methods[data$sample == sample_name]  = as.logical(as.character(toupper(doublet$Summary >= 3)))
    #data$doublet_scores[data$sample == sample_name]  = scrb$doublet_scores
    #data$predicted_doublet_scrub[data$sample == sample_name]  = as.logical(as.character(toupper(scrb$predicted_doublet)))
    
    
  }
  return (data)
}

load_CellLabel <- function(data) {
  data$CellType = ''
  sample_list = unique(data$sample)
  for (i in 1:length(sample_list)){
    #browser()

    sample_name = sample_list[i]
    print(sample_name)

    path = paste0('/home/sujwary/Desktop/scRNA/Output/TestNormalization/Soup_MT_C100/Scran/',sample_name,'/cellIdents.csv')
    cellIdents = read.csv(path,sep = ',',row.names = 1)
    cellIdents$x = paste0(cellIdents$x, ' S',i)
    print(unique(cellIdents$x))
    
    cell_names_all = colnames(data[,data$sample == sample_name])
    cell_names_all = sub("_.*", "", cell_names_all)
    
    cell_names_sample = rownames(cellIdents)
    cell_names_sample = sub("_.*", "", cell_names_sample)
    cellIdents = cellIdents[cell_names_sample %in% cell_names_all,]
    cell_names_sample = cell_names_sample[cell_names_sample %in% cell_names_all]
    
    cellIdents = cellIdents[match(cell_names_sample,cell_names_all)]
    
    data$CellType[data$sample == sample_name] = cellIdents
  }
  data$GeneralCellType = str_match(data$CellType, "(^.+)\\s")[, 2]
  return(data)
}

addMetaData = function(data, metaData){
  #browser()
  #data@meta.data[,colnames(metaData)] = ''
  for (sample in unique(data$sample)){
    print(sample)
    #browser()
    for (colname in colnames(metaData)){
      #print(colname)
      if (colname == 'Dexa or not' & sample == 'GL1080BM'){
        
        #browser()
      }
      #browser()
      values = (as.data.frame(metaData[,colname]))
      values = as.character(values[,1])
      
      metaData_sample_val = values[metaData$Sample == sample]
      metaData_sample_val = metaData_sample_val[[1]]
      if(length(metaData_sample_val) == 0){
        next
      }
      if ( !(colname %in% colnames(data@meta.data)) ){
        data@meta.data[,colname] = ''
      }
      if (!is.na(metaData_sample_val )){
        data@meta.data[,colname][data$sample == sample] = metaData_sample_val
      }
    }
  }
  
  data$LowCount = F
  data$LowCount[data$nCount_RNA < 200] = T
  return (data)
}


