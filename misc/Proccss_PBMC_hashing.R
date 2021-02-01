library(Seurat)
library(Matrix)

library(DropletUtils)
library(ggplot2)
library(SoupX)
source('/home/sujwary/Desktop/scRNA/Code/Plot_func.R')
source('/home/sujwary/Desktop/scRNA/Code/Functions.R')

pbmc.umis <- read.delim("/disk2/Projects/PBMC/PBMC_hashing/GSM2895282_Hashtag-RNA.umi.txt", 
                        header = T, sep = "",row.names = 1)

# Load in the HTO count matrix
pbmc.htos <- read.csv("/disk2/Projects/PBMC/PBMC_hashing/GSM2895283_Hashtag-HTO-count.csv",row.names = 1, header = T)

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis,min.cells = 3, min.features = 1)



pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)



pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
table(pbmc.hashtag$HTO_classification.global)

#
#Save pbmc.hashtag

pbmc_subset <- pbmc.hashtag[,pbmc.hashtag$HTO_classification != 'Negative']

pbmc_subset <- pbmc_subset[,pbmc_subset$HTO_classification != 'total-reads']
pbmc_subset <- pbmc_subset[,pbmc_subset$HTO_classification != 'no-match_total-reads']

pbmc_subset$HTO_classification = gsub("-.*","",pbmc_subset$HTO_classification)
sample_name = unique(pbmc_subset$HTO_classification)[1]

for (sample_name in unique(pbmc_subset$HTO_classification)){
  print(sample_name)
  pbmc_subset_HTO = pbmc_subset[,pbmc_subset$HTO_classification == sample_name]


  pbmc_subset_HTO[["percent.mt"]] <- PercentageFeatureSet(pbmc_subset_HTO, pattern = "^MT-")
  
  dir.create( paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/','QC'), recursive = TRUE)
  pathName <- paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/','QC/',sample_name, '_QC_prefilter.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(pbmc_subset_HTO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0)
  print(plot)
  dev.off()

  data_i_filtered_tmp = pbmc_subset_HTO[,pbmc_subset_HTO$nFeature_RNA< 2000 & pbmc_subset_HTO$percent.mt < 15 ]
  pathName <- paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/','QC/',sample_name, '_QC_maxFeature2000.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data_i_filtered_tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0)
  print(plot)
  dev.off()
  ##
  
  ## Filter
  print(sample_name)
  colSum_list = colSums(pbmc_subset_HTO ) # Needs to be from Matrix library
  keep = colSum_list >= 200
  data_i_filtered = pbmc_subset_HTO[,keep]
  
  
  percent_MT_list = data_i_filtered$percent.mt
  
  pathName <- paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/','QC/',sample_name, '_QC_filter.png')
  png(file=pathName,width=500, height=500)
  plot = VlnPlot(data_i_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3,pt.size = 0)
  print(plot)
  dev.off()
  
  #next

  data_i_filtered_run = NormalizeData(data_i_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
  #data_i_filtered_run = ScranNorm(data_i_filtered)
  data_i_filtered_run = FindVariableFeatures(data_i_filtered_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_run = ScaleData(data_i_filtered_run)
  data_i_filtered_run = RunPCA(data_i_filtered_run,npcs = 30)
  data_i_filtered_run = FindNeighbors(data_i_filtered_run, dims = 1:30)
  data_i_filtered_run = FindClusters(data_i_filtered_run)
  data_i_filtered_run = RunUMAP(data_i_filtered_run, dims = 1:30)
  
  data_matrix_filtered = data_i_filtered_run@assays[["RNA"]]@counts
  data_matrix_raw = pbmc_subset_HTO@assays[["RNA"]]@counts
  cluster_IDs = factor(as.character(Idents(data_i_filtered_run)))
  
  
  sc = SoupChannel(data_matrix_raw, data_matrix_filtered)
  sc = setClusters(sc, cluster_IDs)
  sc = setDR(sc, data_i_filtered_run@reductions[["umap"]]@cell.embeddings)
  
  igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
              "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
  igGenes = igGenes[igGenes %in% sc[["toc"]]@Dimnames[[1]]] # usetoEst will break if not all genes are present
  HBGenes = c('HBB','HBA2')
  HBGenes = HBGenes[HBGenes %in% sc[["toc"]]@Dimnames[[1]]] 
  
  useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = list(IG = igGenes))
  
  
  sc = calculateContaminationFraction(sc, list(IG = igGenes), useToEst = useToEst)
  
  folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PreSoup_igGenes','','.png')
  png(file=pathName,width=1000, height=1000)
  print( plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst))
  dev.off()
  
  folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  
  
  out = adjustCounts(sc)
  print('Plot')
  pathName = paste0(folder,sample_name,'_PreSoup_Umap','_SoupAdjust','.png')
  png(file=pathName,width=1000, height=1000)
  print(plotChangeMap(sc, out, geneSet = HBGenes))
  dev.off()
  
  
  # Put soup data back into filtered run
  data_i_filtered_soup = data_i_filtered
  
  data_i_filtered_soup@assays[["RNA"]]@counts = out
  
  data_i_filtered_soup = data_i_filtered_soup[, data_i_filtered_soup$percent.mt < 15]
  
  data_i_filtered_soup_run = NormalizeData(data_i_filtered_soup, normalization.method = "LogNormalize", scale.factor = 10000)
  #data_i_filtered_run = ScranNorm(data_i_filtered)
  data_i_filtered_soup_run = FindVariableFeatures(data_i_filtered_soup_run, selection.method = "vst", nfeatures = 2000)
  data_i_filtered_soup_run = ScaleData(data_i_filtered_soup_run)
  data_i_filtered_soup_run = RunPCA(data_i_filtered_soup_run,npcs = 30)
  data_i_filtered_soup_run = FindNeighbors(data_i_filtered_soup_run, dims = 1:30)
  data_i_filtered_soup_run = FindClusters(data_i_filtered_soup_run)
  data_i_filtered_soup_run = RunUMAP(data_i_filtered_soup_run, dims = 1:30)
  
  
  
  folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/')
  dir.create(folder, recursive = T)
  pathName = paste0(folder,sample_name,'_Post_Soup_Umap','','.png')
  png(file=pathName,width=1000, height=1000)
  print(  DimPlot(data_i_filtered_soup_run,pt.size = 0.5, reduction = "umap",label = FALSE))
  dev.off()
  
  path = paste0(folder,'/',sample_name,'.Robj')
  save(data_i_filtered_soup_run,file= path)
  
  write(colnames(data_i_filtered_soup@assays[["RNA"]]@counts ), file = paste0(folder,sample_name,'_colnames.txt'))
  write(rownames(data_i_filtered_soup@assays[["RNA"]]@counts ), file = paste0(folder,sample_name,'_rownames.txt'))
  writeMM(data_i_filtered_soup@assays[["RNA"]]@counts , file = paste0(folder,sample_name,'_matrix.txt'))
  
  next
  gene_list = c('CD3D', 'CD3G', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'SELL', 'CD14', 'FCGR3A', 'NKG7', 'MS4A1', 'IGKC', 'IGHM', 'CD19', 'MZB1', 'CD34', 'CDK6',
                'FCER1A','FUT4', 'ELANE', 'MPO', 'HBA2', 'HBB', 'LYZ', 'TNFRSF17')
  
  for (j in 1:length(gene_list)){
    gene = gene_list[j]
    print(gene)
    #browser()
    folder = paste0('/disk2/Projects/EloRD/Output/Soup_MT_C100/',sample_name,'/','Featureplots/')
    dir.create(folder,recursive = T)
    plot = FeaturePlotFix(data_i_filtered_soup_run, feature = gene, folder =folder,
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

## Empty drops
for (sample_name in unique(pbmc_subset$HTO_classification)){
  
  file = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/', 'emptyDrops','.csv')
  if(file.exists( file)){
    print('File already exists')
    next
  }
  #Scrublet_threshold = sampleParam$Scrublet_threshold[sampleParam['Sample'] == sample_name]
  #print(Scrublet_threshold)
  
  print(sample_name)
  pbmc_subset_HTO = pbmc_subset[,pbmc_subset$HTO_classification == sample_name]
  
  
  pbmc_subset_HTO[["percent.mt"]] <- PercentageFeatureSet(pbmc_subset_HTO, pattern = "^MT-")
  
  ## Empty drops
  counts= pbmc_subset_HTO@assays[["RNA"]]@counts
  br.out = barcodeRanks(counts,lower=100)
  e.out = emptyDrops(counts,lower=100) # Outputs NAs for cell with lower than 100
  
  
  br = data.frame(br.out@rownames,br.out$rank, br.out$total)
  e = data.frame(e.out@rownames, e.out$FDR,e.out$LogProb)
  colnames(br) = c('rownames','rank','total')
  colnames(e) = c('rownames','FDR','LogProb')
  
  br_e <- merge(br,e,by="rownames")
  br_e$is_cell <- br_e$FDR <= 0.01
  br_e$ncount = pbmc_subset_HTO$nCount_RNA
  
  #ggplot(br_e) + geom_point(aes(rank,total,color=is_cell)) + scale_y_log10() + scale_x_log10()
  
  
  br_e_sorted = br_e[match(as.character(rownames(data_i_filtered@meta.data)), as.character(br_e$rownames)   ),]
  
  folder = paste0('/home/sujwary/Desktop/scRNA/Output/EmptyCells/',sample_name,'/')
  dir.create(folder, recursive =T)
  write.csv(br_e, file = paste0(folder,'emptyDrops','.csv'),row.names = FALSE)
}
###############

pbmc_subset <- pbmc.hashtag[,pbmc.hashtag$HTO_classification.global != 'Negative']
pbmc_subset <- pbmc_subset[,pbmc_subset$HTO_classification != 'no-match_total-reads']
pbmc_subset <- pbmc_subset[,pbmc_subset$HTO_classification != 'total-reads']


# Select the top 1000 most variable features
pbmc_subset <- FindVariableFeatures(pbmc_subset, selection.method = "vst")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc_subset <- ScaleData(pbmc_subset, features = VariableFeatures(pbmc_subset))

# Run PCA
pbmc_subset <- RunPCA(pbmc_subset, features = VariableFeatures(pbmc_subset))
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc_subset <- FindNeighbors(pbmc_subset, reduction = "pca", dims = 1:10)
pbmc_subset <- FindClusters(pbmc_subset, resolution = 0.6, verbose = FALSE)
pbmc_subset <- RunTSNE(pbmc_subset, reduction = "pca", dims = 1:10,check_duplicates = F)

# Projecting singlet identities on TSNE visualization
DimPlot(pbmc_subset, group.by = "HTO_classification")

