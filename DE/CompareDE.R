
#######################
## DE
#######################
# Detection rate
#data = as.data.frame(data_harmony_run_label@assays[["RNA"]]@counts)
#DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
#data_harmony_run_label$DR = DR
source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')

data_harmony_run_label = ScaleData(data_harmony_run_label, features = rownames(data_harmony_run_label))


celltype_list = unique(as.character(Idents(data_harmony_run_label)))
celltype_list

celltype_num = sort(table(as.character(Idents(data_harmony_run_label))))

celltype_list = c('cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')

celltype_list = c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell')


celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')



celltype_list = c('Cytotoxic NK')

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')

celltype_list = c('dMono','Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype_list = c('GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype = celltype_list[1]
celltype_list = c('Inhibitory NK')


#celltype_list = c('Plasma Cell')
#celltype_list=  c('CD56Br NK','Th17','Inhibitory NK','TEMRA','DC')
for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0('baseline Yes ',celltype)
  ident2 = paste0('baseline No ',celltype)
  
  #ident1 = paste0('CD16+ Mono')
  #ident2 = paste0('CD14+ Mono')
  
  DE_input = data_harmony_run_label
  DE_input = renameCells(DE_input,idents = c('cDC1','cDC2'),newident = 'DC')
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  #DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
  #                          DE_input$Best_Overall_Response, ' ', Idents(DE_input))
  DE_input$DE_ident = paste0(DE_input$Treatment,' ',DE_input$'Dexa or not', ' ', Idents(DE_input))
  #DE_input$DE_ident = paste0(Idents(DE_input))
  #DE_input = DE_input[,DE_input$Treatment =='baseline']
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
  ncol(DE_input)
  unique(DE_input$DE_ident)
  
  if (ncol(DE_input) < 100){
    next
  }
  
  #DE_input_sce = as.SingleCellExperiment(DE_input)
  
  data = as.data.frame(DE_input@assays[["RNA"]]@counts)
  # detection rate:fraction of genes expressed in a cell
  DR = as.numeric(apply(data,2, function(x) sum(x > 0)/nrow(data)))
  DE_input$DR = DR 
  
  kit = factor(DE_input$kit)
  ident = factor(as.character(Idents(DE_input)))
  DE_input$ident = ident
  DE_ident = factor(DE_input$DE_ident)
  Patient = factor(DE_input$`Patient Number`)
  DE_input$Patient = Patient
  Treatment = factor(DE_input$Treatment)
  
  if (length(unique(DE_ident)) == 1){
    next
  }
  library(robustbase)
  library(DESeq2)
  
  if (length(unique(kit)) == 1){
    formula = ~  DE_ident + DR + Patient
  }else{
    formula = ~  DE_ident + DR + kit + Patient
  }
  #formula = ~  DE_ident + DR + kit
  design <- model.matrix(formula)
  colnames(design) <- gsub("DE_ident", "", colnames(design))
  colnames(design)
  is.fullrank(design)
  design_fr = fullRank(design)
  is.fullrank(design_fr)
  colnames_old = colnames(design)
  colnames_fr = colnames(design_fr)
  
  colnames_old[!(colnames_old %in% colnames_fr)]
  
  keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                       min.count = 1,min.total.count=10, 
                       large.n = 10,min.prop = 0.1)
  subfolder = paste0(ident1,' Vs ',ident2)
  
  result_DESeq2 = runDESeq2(DE_input,design_fr,contrast, keep, 
                            folder_output = filepath_cluster, subfolder = subfolder)
  
  dds = result_DESeq2[[1]]
  res_DESeq2 = result_DESeq2[[2]]
  #plotDispEsts(dds)
  
  res_DESeq2 = res_DESeq2[order(res_DESeq2$log2FoldChange),]
  res_DESeq2 = res_DESeq2[res_DESeq2$pvalue < 0.05,]
  
  res_DESeq2 = res_DESeq2[!grepl("MT-", res_DESeq2$gene),]
  res_DESeq2 = res_DESeq2[!grepl(" ?RP\\w+ ?", res_DESeq2$gene),]
  
  folder = paste0(filepath_cluster,'/DE/DESeq2/',subfolder,'/')
  dir.create(folder,recursive = T)
  pathName <- paste0(folder,paste0('HeatMap','.png'))
  #png(file=pathName, height=12, width=20)
  #png(file=pathName)
  len = nrow(res_DESeq2)
  if (nrow(res_DESeq2) > 200){
    res_DESeq2 = res_DESeq2[c(1:100,(len - 100):len),]
    res_DESeq2 = res_DESeq2[rowSums(is.na(res_DESeq2)) != ncol(res_DESeq2),]
    res_DESeq2 = unique(res_DESeq2)
  }
  
  
  png(file=pathName,width=2000, height=40*nrow(res_DESeq2),res = 100)
  
  plot = DoHeatmap(DE_input, features = res_DESeq2$gene, group.by = 'DE_ident')
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=24 ),
    axis.title.y = element_text(color="black", size=24),
    axis.text= element_text(color="black", size=24),
    legend.text=element_text(size=24),
    legend.title=element_text(size=24),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
}
########################



#######################
sum(res_DESeq2$padj < 0.05, na.rm = T)
sum(res_DESeq2$padj > 0.05, na.rm = T)
nrow(res_DESeq2)

res_DESeq2_sig = res_DESeq2$gene[res_DESeq2$padj < 0.05 & !is.na(res_DESeq2$padj)]



res <- results(dds, 
               name = names[2],
               alpha = 0.05,  pAdjustMethod	 ='BH')
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_sig = res_tbl[res_tbl$padj < 0.05 & !(is.na(res_tbl$padj)),]
res_tbl_sig


Idents(DE_input) = DE_input$DE_ident
Features = FindMarkers(DE_input, ident.1 = ident1, ident.2 = ident2
                       ,min.pct = 0.1, logfc.threshold = 0.1, only.pos = FALSE, slot = "counts")

Features_sig = Features[Features$p_val_adj< 0.05,]

intersect(res_tbl_sig$gene, rownames(Features_sig))
#  Identifies differentially expressed genes between two groups of cells 
#  using a hurdle model tailored to scRNA-seq data. 
# Utilizes the MAST package to run the DE testing.

# MAST
markers = FindMarkers(data_harmony_run_label,  
                      ident.1 = "CD8+ T Cell", ident.2 = "T Cell", 
                      latent.vars = 'kit', test.use = "MAST")

##
# edgeR
data_harmony_run_label_sce = as.SingleCellExperiment(data_harmony_run_label)
## Convert to DGEList, calculate logCPMs
dge <- scran::convertTo(data_harmony_run_label_sce, type = "edgeR")
plotMDS(dge)

design <- model.matrix(~kit+DR)

dge <- estimateDisp(dge, design, robust=TRUE)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)
topTags(lrt)

##
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)



# voom limma
dge <- DGEList(data_harmony_run_label@assays[["RNA"]]@counts, 
               group = data_harmony_run_label$kit)
dge <- calcNormFactors(dge)

design <- model.matrix(~ kit + DR)

vm <- voom(dge, design = design, plot = TRUE)
fit <- lmFit(vm, design = design)
fit <- eBayes(fit)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

hist(tt$P.Value, 50)
hist(tt$adj.P.Val, 50)
#limma::plotMDS(dge, col = as.numeric(as.factor(data_harmony_run_label$kit)), pch = 19)
plotMD(fit)

# Cluster comparisions using design matrix made from cell idetents
design <- model.matrix(~ 0 + ident, data = colData(data_harmony_run_label_sce))
colnames(design) <- gsub("ident", "", colnames(design))
colnames(design)

y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- limma::lmFit(y, design)

## Perform pairwise comparisons
nclust <- length(unique(data_harmony_run_label_sce$ident))
all.results <- all.pairs <- list()
counter <- 1

for (i in seq_len(nclust)) {
  for (j in seq_len(i - 1L)) {
    con <- integer(ncol(design))
    con[i] <- 1
    con[j] <- -1
    fit2 <- limma::contrasts.fit(fit, con)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    
    res <- limma::topTable(fit2, number = Inf, sort.by = "none")
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(i, j)]
    counter <- counter + 1L
    
    ## Also filling the reverse comparison.
    res$logFC <- -res$logFC
    all.results[[counter]] <- res
    all.pairs[[counter]] <- colnames(design)[c(j, i)]
    counter <- counter + 1L
  }
}

## Combine results across all pairwise tests
all.pairs <- do.call(rbind, all.pairs)
combined <- scran::combineMarkers(all.results, all.pairs, 
                                  pval.field = "P.Value",
                                  pval.type = "any")
head(combined[["cluster1"]])

markers <- scran:::findMarkers(data_harmony_run_label_sce, design=design)
