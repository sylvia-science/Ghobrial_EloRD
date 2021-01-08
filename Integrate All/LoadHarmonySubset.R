LoadSubset = function(data_merge_run_label,sampleParam_combine, folder){
  
  browser()
  celltype = 'Mono_DC'
  resolution_val_subset = 1.6
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  
  Ident_main = colnames(data_merge_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_merge_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_merge_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_merge_run_label) = newIdents2
  ##
  
  celltype = 'NK'
  resolution_val_subset = 3
  cluster_IDs_subset =sampleParam_combine$Cluster_IDs[sampleParam_combine['Sample'] == paste0(sample_type,'_',celltype)]
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  
  Ident_main = colnames(data_merge_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_merge_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_merge_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_merge_run_label) = newIdents2
  
  
  #
  
  celltype = 'T Cell'
  resolution_val_subset = 3.5
  
  
  folder_subcluster = paste0(folder, 'Subcluster/',celltype,'/')
  path = paste0(folder_subcluster,'data_run','.Robj')
  data_run_subset = loadRData(path)
  tmp = data_run_subset@meta.data[paste0('RNA_snn_res.', resolution_val_subset)]
  tmp = tmp[,1]
  tmp  =factor(tmp,levels = as.character(0:(length(unique(tmp))-1)))
  Idents(data_run_subset) = tmp
  
  data_run_subset_label =label_cells(data_run_subset,cluster_IDs_subset)
  
  
  Ident_main = colnames(data_merge_run_label)
  Ident_main = Ident_main[Ident_main %in% colnames(data_run_subset_label)]
  
  Ident_subset = colnames(data_run_subset_label)
  Ident_subset_match = match(Ident_subset, Ident_main)
  cell_subset = Idents(data_run_subset_label)
  
  newIdents = as.character(Idents(data_merge_run_label))
  newIdents2= newIdents
  newIdents2[colnames(data_merge_run_label) %in% colnames(data_run_subset_label) ] = as.character(cell_subset[Ident_subset_match])
  Idents(data_merge_run_label) = newIdents2
  
  
  return(data_merge_run_label)
  
}
  