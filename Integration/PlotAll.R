plotAll = function(data,folder,sample_name,sampleParam,label_TF,
                   cell_features,
                   plot_PCA= F,
                   DE_perm_TF = FALSE, # Do DE between each cluster and each individual cluster
                   DE_perm_clusters = NA,
                   markersTF = T, # Do DE between each cluster and all other clusters combined.
                   groupBy = NA, # Variables that will be used to make umap plots where the variables are split in the plot
                   splitBy = NA, # Variables to plot split in the umap
                   featurePlot_list = NA,
                   PCA_dim = NA, # If NA, will gram PCA_dim from sampleParam
                   resolution_val = NA, # If NA, will get resolution_val from sampleParam
                   pt.size = 0.7,
                   str = '' 
){
  require(gtools)
  #browser()
  print(folder)
  print(paste0('label_TF: ', label_TF))
  
  data$percentMT = data$percent.mt
  #s.genes <- cc.genes$s.genes
  #g2m.genes <- cc.genes$g2m.genes
  #data = CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE,group.singletons = FALSE)
  
  if (is.na(PCA_dim)){
    Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
    PCA_dim = sampleParam$PCA_dim[sampleParam['Sample'] == sample_name]
  }
  if (is.na(resolution_val)){
    Patient_num  = sampleParam$`Patient Number`[sampleParam['Sample'] == sample_name]
    resolution_val = sampleParam$resolution_val[sampleParam['Sample'] == sample_name]
  }
  
  filepath_cluster = paste0( folder, 'Cluster/', 'PCA',PCA_dim,'/res',resolution_val,'/' )
  print('filepath_cluster')
  print(filepath_cluster)
  dir.create( filepath_cluster, recursive = TRUE)
  
  # Score for cell cycle
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  #browser()
  
  if (label_TF){
    file_str = '_label'
  }else
  {
    file_str = ''
  }
  file_str = paste0(file_str, str)
  print(paste0('file_str: ', file_str))
  ##################
  ## QC
  ##################
  
  OrigIdents = Idents(data)
  Idents(data) = ''
  dir.create( paste0(folder,'QC'), recursive = TRUE)
  pathName <- paste0(folder,'QC/', 'QC.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0)
  print(plot)
  dev.off()
  Idents(data) = OrigIdents
  
  
  ##################
  ## PCA
  ##################
  #browser()
  reduction = 'pca'
  reduction = 'harmony'
  if (plot_PCA){
    visualize_PCA(data,folder,PCA_dim,reduction)
    dir.create( paste0(folder,reduction), recursive = TRUE)
    pathName <- paste0(folder,reduction,'/elbow_',PCA_dim,'.png')
    png(file=pathName,width=600, height=350)
    print(ElbowPlot(data,ndims = PCA_dim))
    dev.off()
  }
  #############
  ## TEST: REFIND VAR FEATURE
  #############
  
  #browser()
  nfeatures_val = 2000 
  #data_RNA = data
  #DefaultAssay(object = data_RNA) <- "RNA"
  #data_RNA@active.assay = 'RNA'
  #data_post_varFeatures = FindVariableFeatures(data_RNA, selection.method = "vst", nfeatures = nfeatures_val)
  
  #browser()
  
    #data = FindNeighbors(data, dims = 1:PCA_dim)
    #data = FindClusters(data, resolution = resolution_val)
  
  #browser()
  
  # Name cells
  if (label_TF){
    #browser()
    cluster_IDs <- sampleParam$Cluster_IDs[sampleParam['Sample'] == sample_name]
    data = label_cells(data,cluster_IDs)
    # Percentage of cell type
  }
  
  ##############################
  ## Get Cluster Stats
  ##############################
  
  #cluster_num = clusterStats(data)
  #print(cluster_num)
  
  #filepath_stats = paste0(filepath_cluster, 'Stats/')
  #dir.create( filepath_stats, recursive = TRUE)
  #write.csv(cluster_num, file = paste0(filepath_stats,'clusterStats',file_str,'.csv'),row.names = FALSE)
  
  
  for (var in unique(data$split_var)){
    print(var)
    #data_subset = data[,data$split_var == var]
    #cluster_num = clusterStats(data_subset)
    #print(cluster_num)
    
    #filepath_stats = paste0(filepath_cluster, 'Stats/')
    #dir.create( filepath_stats, recursive = TRUE)
    #write.csv(cluster_num, file = paste0(filepath_stats,'clusterStats_',var,file_str,'.csv'),row.names = FALSE)
    
  }
  
  
  ########################
  ## Visualize clustering
  ########################
  
  #browser()
  pathName <- paste0(filepath_cluster,
                     paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
  png(file=pathName,width=2500, height=1500, res = 100)
  plot = DimPlot(data,pt.size = pt.size, reduction = "umap",
                 label = TRUE,label.size = 14)
  plot = plot + theme(
    legend.title = element_text( size = 24),
    legend.text = element_text( size = 24))
  plot = plot +theme(axis.text=element_text(size=24),
                     axis.title=element_text(size=24,face="bold"))
  print(plot)
  dev.off()
  
  pathName <- paste0(filepath_cluster,
                     paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_No_label',file_str,'.png'))
  png(file=pathName,width=2500, height=1500, res = 100)
  plot = DimPlot(data,pt.size = pt.size, reduction = "umap",
                 label = F,label.size = 14)
  plot = plot + theme(
    legend.title = element_text( size = 24),
    legend.text = element_text( size = 24))
  plot = plot +theme(axis.text=element_text(size=24),
               axis.title=element_text(size=24,face="bold"))
  
  print(plot)
  dev.off()
  
  
  #pathName <- paste0(filepath_cluster,'ClusterMetrics','.png')
  #png(file=pathName,width=600, height=350)
  #print(FeaturePlot(data,pt.size = 0.5, features = c("percent.mt")))
  #dev.off()

  
  pathName <- paste0(filepath_cluster,'nCount_RNA','.png')
  png(file=pathName,width=600, height=600)
  print(FeaturePlot(data,pt.size = pt.size, features = c("nCount_RNA")))
  dev.off()
  
  
  #pathName <- paste0(filepath_cluster,'S.Score','.png')
  #png(file=pathName,width=600, height=600)
  #print(FeaturePlot(data,pt.size = 0.5, features = c("S.Score")))
  #dev.off()

 
  
  
  #browser()
  
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap','_PCA',PCA_dim,resolution_val,'_splitAll', '','.png'))  
  png(file=pathName,width=2600, height=500,res = 100)
  print(DimPlot(data, label=T, repel=F, reduction = "umap", split.by = "split_var"))
  dev.off()
  
  #browser()
  if (!is.na(featurePlot_list)){
    for (feature in featurePlot_list){
      pathName <- paste0(filepath_cluster,feature,'.png')
      png(file=pathName,width=600, height=600)
      print(FeaturePlot(data,pt.size = pt.size, features = feature))
      dev.off()
    }
  }
  if (!is.na(groupBy)){
    for (group in groupBy){
      print(group)
      pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy',group,'.png'))
      png(file=pathName,width=1000, height=1000)
      plot = DimPlot(data,pt.size = pt.size, reduction = "umap",
                     label = F,group.by  = group)
      
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20)
      )
      print(plot)
      
      dev.off()
      
     
      
    }
  }
  
  if (!is.na(splitBy)){
    for (group in splitBy){

      pathName <- paste0(filepath_cluster,paste0('ClusterUmap','_PCA',PCA_dim,'_res',resolution_val,'_split', group,'.png'))  
      png(file=pathName,width=3000, height=1000)
      plot = DimPlot(data, pt.size = pt.size, 
                     label=T, repel=F, reduction = "umap", split.by = group)
      plot = plot + theme(
        axis.title.x = element_text(color="black", size=24 ),
        axis.title.y = element_text(color="black", size=24),
        axis.text= element_text(color="black", size=24),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24),
        text = element_text(size = 20)
        
      )
      print(plot)
      dev.off()
      
    }
  }
  #browser()
  
  # Plot larger sample plot
  group = 'sample'
  pathName <- paste0(filepath_cluster,paste0('ClusterUmap', '_PCA',PCA_dim,'_res',resolution_val,'_GroupBy',group,'_wide.png'))
  png(file=pathName,width=3000, height=1000)
  plot = DimPlot(data,pt.size = pt.size, reduction = "umap",
                 label = F,group.by  = group)
  
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20)
  )
  print(plot)
  dev.off()
  ##################################
  
  ##############################
  ## Find Cluster Biomarkers
  ##############################
  cell_features_label = cell_features[cell_features$Use_Label == 1,]
  
  if (markersTF == TRUE){
    # find markers for every cluster compared to all remaining cells, report only the positive ones
    num_markers = 10
    markers = FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers  %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
    # Plotting the top 10 markers for each cluster.
    top10 = markers %>% group_by(cluster) %>% top_n(n = num_markers, wt = avg_logFC)
    all_markers =  markers %>% group_by(cluster)
    #browser()
    all_markers = cellMarkers(all_markers,cell_features_label)
    write.csv(all_markers, file = paste0(filepath_cluster,'Features',file_str,'.csv'),row.names=FALSE)
    
    pathName <- paste0(filepath_cluster,paste0('HeatMap', '_PCA',PCA_dim,'_res',resolution_val,file_str,'.png'))
    png(file=pathName,width=1000, height=1500)
    plot = DoHeatmap(data, features = top10$gene)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=18 ),
      axis.title.y = element_text(color="black", size=18),
      axis.text= element_text(color="black", size=18),
      legend.text=element_text(size=18),
      legend.title=element_text(size=18),
      text = element_text(size = 18))
    print(plot)
    dev.off()
    #browser()

  } 
  
  # Get DE between cluster permutations
  if (DE_perm_TF){
    #browser()
    filepath_mvn = paste0( filepath_cluster, 'DE/nVsm/')
    print(filepath_mvn)
    dir.create( filepath_mvn, recursive = TRUE)
    if (!(is.na(DE_perm_clusters))){
      data_input = data[,Idents(data) %in% DE_perm_clusters]
    }
    #browser()
    Features_mvn_df = diffExpPermute(data_input,cell_features_label)
    
    level_list = unique(levels(data_input))
    for (level in 1:length(level_list) ){
      ident1 = level_list[level]
      
      if (ident1 == '?'){
        ident1 = 'NA'
      }
      
      ident1 = gsub("/", "", ident1, fixed = T)

      df_output = Features_mvn_df[[level]]
      df_output = df_output[c("p_val","avg_logFC", "pct.1", "pct.2", "p_val_adj","ident_1","ident_2","gene", "Cell")]
      filepath = paste0(filepath_mvn
                        ,'Features_',ident1,'Vsn',file_str
                        ,'.csv')
      print(filepath)
      #browser()
      df_output = cellMarkers(df_output,cell_features_label)
      
      write.csv(df_output, file = filepath,row.names = FALSE)
    }
  }
  
}


