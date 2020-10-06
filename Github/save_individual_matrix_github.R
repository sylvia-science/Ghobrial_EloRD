library(Matrix)
library(readxl)
library(Seurat)

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
metaData = metaData[metaData$Run== 1,]

sample_list = metaData$Sample
for (i in 1:length(sample_list)){
  sample_name = sample_list[i]
  print(sample_name)
  #sample_name = 'GL1080BM' #metaData$Sample[4]
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_raw_feature_bc_matrix.h5",sep = "")
  
  #filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample_name,"_filtered_feature_bc_matrix.h5",sep = "")
  path = paste0('/home/sujwary/Desktop/scRNA/Output/Doublets/',sample_name,'_scrublet0.3.Robj')
  
  
  data_i_scrublet = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  #data_i_scrublet = CreateSeuratObject(counts = data_i_scrublet, project = "BM", min.cells = 3,min.features = 1)
  #colSum_list = colSums(data_i_scrublet)
  #keep = colSum_list >= 100
  
  #data_i_scrublet = data_i_scrublet[,keep]
  
  #data_matrix = data_i_scrublet@assays[["RNA"]]@counts
  data_matrix = data_i_scrublet
  
  write.csv(colnames(data_matrix), file = paste0('/disk2/Projects/EloRD/Data/RawMatrix/',sample_name,'_out_cell_barcodes_filt.csv'),
                                                 quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #write(rownames(data_matrix), file = paste0('/disk2/Projects/EloRD/Data/RawMatrix',sample_name,'_rownames.txt'))
  #writeMM(data_matrix, file = paste0('/disk2/Projects/EloRD/Data/RawMatrix/',sample_name,'_matrix.txt'))
  
}
  
