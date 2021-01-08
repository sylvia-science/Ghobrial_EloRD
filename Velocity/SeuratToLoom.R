dir.create(paste0(filepath_cluster,'Data/'))
write.csv(as.character(data_run_subset$cell_sample), file = paste0(filepath_cluster,'Data/',"cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(data_run_subset, reduction = "umap"), file = paste0(filepath_cluster,'Data/',"cell_embeddings.csv"))
write.csv(data_run_subset@meta.data$seurat_clusters, file =paste0(filepath_cluster,'Data/', "clusters.csv"))

##

data_input = data_run_subset

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


library(Seurat)
library(velocyto.R)
library(SeuratWrappers)

bm <- RunVelocity(object = data_input, deltaT = 1, kCells = 25, fit.quantile = 0.02,spliced = 'RNA',unspliced  = 'RNA')
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
                               slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)
