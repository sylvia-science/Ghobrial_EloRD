


library(reticulate)
use_condaenv("r-velocity", required = T)

reticulate::py_config()
conda_list()


scv <- import("scvelo")
scv$logging$print_version()

#####################
adata <- scv$datasets$pancreas()
adata


#scv$pl$scatter(adata, legend_loc='lower left', size=60)

## get embedding
emb <- adata$obsm['X_umap']
clusters <- adata$obs$clusters
rownames(emb) <- names(clusters) <- adata$obs_names$values

## get clusters, convert to colors
col <- rainbow(length(levels(clusters)), s=0.8, v=0.8)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)

## simple plot
plot(emb, col=cell.cols, pch=16,
     xlab='UMAP X', ylab='UMAP Y')
legend(x=-13, y=0, 
       legend=levels(clusters),
       col=col, 
       pch=16)

## run scvelo dynamic model
scv$pp$filter_genes(adata) ## filter
scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model

## takes awhile, so uncomment to save
#adata$write('data/pancreas.h5ad', compression='gzip')
#adata = scv$read('data/pancreas.h5ad')

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
scv$pl$velocity_embedding_stream(adata, basis='umap')
## scv$pl$velocity_embedding_stream(adata, basis='pca') ## other embedding

## top dynamic genes
topgenes <- adata$var["fit_likelihood"]
topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)

scv$pl$scatter(adata, basis=names(topgenes_vals)[1:5], ncols=5, frameon=FALSE)


pcs <- adata$obsm['X_pca']
rownames(pcs) <- cells
##
## Own data
##
data_run_subset

data_run_subset_anndata <- Convert(data_run_subset, to = "anndata")


ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(Idents(data_run_subset), Idents(data_run_subset))
rownames(dfobs) <- colnames(data_run_subset)
dfvar <- data_run_subset@meta.data
rownames(dfvar) <- rownames(data_run_subset)

pcs = data_run_subset@reductions[["pca"]]@cell.embeddings
umap = data_run_subset@reductions[["umap"]]@cell.embeddings

adata_input <- ad$AnnData(
  X=t(data_run_subset@assays[["RNA"]]@counts),
  obs=dfobs,
  obsm=list('X_umap'=umap, 'X_pca'=pcs[,1:30]) 
)
adata_input$obs$clusters = Idents(data_run_subset)

#scv$pl$scatter(adata_input, legend_loc='lower left', size=60)

## get embedding
emb <- umap
clusters <- Idents(data_run_subset)
rownames(emb) <- as.character(Idents(data_run_subset))

## get clusters, convert to colors
col <- rainbow(length(levels(clusters)), s=0.8, v=0.8)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)

## simple plot
plot(emb, col=cell.cols, pch=16,
     xlab='UMAP X', ylab='UMAP Y')
legend(x=-13, y=0, 
       legend=levels(clusters),
       col=col, 
       pch=16)

## run scvelo dynamic model
scv$pp$filter_genes(adata_input) ## filter
scv$pp$moments(adata_input) ## normalize and compute moments
scv$tl$recover_dynamics(adata_input) ## model



scv$pp$filter_and_normalize(adata_input, enforce = T)
scv$pp$moments(adata_input)
scv$tl$velocity(adata_input, mode = "stochastic")
scv$tl$velocity_graph(adata_input)
scv$pl$velocity_embedding(adata_input, basis='umap')



## takes awhile, so uncomment to save
#adata$write('data/pancreas.h5ad', compression='gzip')
#adata = scv$read('data/pancreas.h5ad')

## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
scv$tl$velocity_graph(adata)
scv$pl$velocity_embedding_stream(adata, basis='umap')
## scv$pl$velocity_embedding_stream(adata, basis='pca') ## other embedding

## top dynamic genes
topgenes <- adata$var["fit_likelihood"]
topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)

scv$pl$scatter(adata, basis=names(topgenes_vals)[1:5], ncols=5, frameon=FALSE)


pcs <- adata$obsm['X_pca']
rownames(pcs) <- cells

