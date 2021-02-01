data_input = data_run_subset_label_remove
data_input = ScaleData(data_input, features = rownames(data_input))
data_input = renameCells(data_input,idents = c('cDC1','cDC2'),newident = 'DC')
data_input = renameCells(data_input,idents = c('CD14+CD16+ Mono'),newident = 'CD16+ Mono')

data_input$Treatment[data_input$Treatment == 'baseline'] = 'Baseline'

plot = DimPlot(data_input,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

library(fgsea)
pathways.hallmark <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/h.all.v7.2.symbols.gmt')
pathways.c2 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c2.all.v7.2.symbols.gmt')
pathways.c5 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c5.all.v7.2.symbols.gmt')
pathways.c7 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c7.all.v7.2.symbols.gmt')
pathways.c8 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c8.all.v7.2.symbols.gmt')
pathways = c(pathways.hallmark,pathways.c2,pathways.c5, pathways.c7,pathways.c8)


font_size= 32

label_size = 12

gene_height = 80

celltype_list = sort(unique(as.character(Idents(data_input))))
for (celltype in celltype_list){
  print(celltype)
}


DE_type = 'DESeq2'
DE_type = 'EdgeR'

time1 = 'Baseline'
time2 = 'NBM'

base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/', time1,' Vs ',time2,'/')
base = paste0(filepath_cluster,'/DE/',DE_type,'/')

#base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/C9D1 Dexa Vs C9D1 No Dexa/')



celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','Plasma Cell')


celltype_list=  c('NFkB-high')
celltype = celltype_list[1]

boxplot_TF = F

for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0(time1,' ',celltype)
  ident2 = paste0(time2,' ',celltype)
  
  #ident1 = paste0('C9D1 Yes ',celltype)
  #ident2 = paste0('C9D1 No ',celltype)
  
  subfolder = paste0(ident1,' Vs ',ident2)
  
  
  
  DE_input = data_input
  #DE_input = renameCells(DE_input,idents = c('cDC1','cDC2'),newident = 'DC')
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  #DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
  #                          DE_input$Best_Overall_Response, ' ', Idents(DE_input))
  DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
  #DE_input$DE_ident = paste0(DE_input$Treatment,' ',DE_input$'Dexa or not', ' ', Idents(DE_input))
  
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
 
  DE_input$DE_ident = factor(as.character(DE_input$DE_ident), levels =c(ident1,ident2))
  
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'.csv')
  
  if(!file.exists(filename)){
    next
  }
    
  res = read.csv(filename)
  res$gene = res$X
  names(res)[names(res) == "logFC"] <- "log2FoldChange"
  names(res)[names(res) == "FDR"] <- "padj"
  names(res)[names(res) == "PValue"] <- "pvalue"
  
  res = res[order(res$log2FoldChange),]
  res = res[res$padj < 0.05,]
  
  
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'_clean','.csv')
  
  write.csv(res, file = filename,row.names=TRUE)
  
  #next
  
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl("^RP[SL]", res$gene),]
  
  res = res[rowSums(is.na(res)) != ncol(res),]
  
  
  if ('TSC22D3' %in% res$gene){
    print('Found')
    print(celltype)
    print(res[res$gene == 'TSC22D3',])
    #browser()
  }
  #next
  #next
  
  
  ranks <- log10(res$pvalue)/sign(res$log2FoldChange)
  names(ranks) <- res$gene
  fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
  fgseaRes  = fgseaRes[order(fgseaRes$padj, fgseaRes$pval),]
  fgseaRes$leadingEdge =  sapply( fgseaRes$leadingEdge , paste0, collapse=",")
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,'_','fgseaRes_',subfolder,'.csv')
  write.csv(fgseaRes, file = filename,row.names=F)
  
  
  max_gene = 50
  res$gene = as.character(res$gene)
  res_pos = res[res$log2FoldChange > 0,]
  res_neg = res[res$log2FoldChange < 0,]
  if (nrow(res_pos) > max_gene){
  
    res_pos = res_pos[1:max_gene,]
    res_pos = res_pos[rowSums(is.na(res_pos)) != ncol(res_pos),]
    
  }
  if (nrow(res_neg) > max_gene){
    
    res_neg = res_neg[( nrow(res_neg) - max_gene):nrow(res_neg),]
    res_neg = res_neg[rowSums(is.na(res_neg)) != ncol(res_neg),]
    
  }
  res = res[res$gene %in% c(res_pos$gene, res_neg$gene),]
  res = unique(res)
  res = res[rowSums(is.na(res)) != ncol(res),]
  
  res_pos = res_pos[order(-res_pos$log2FoldChange),]
  
  print('pos')
  print(res$gene[res$log2FoldChange > 0])
  
  print('neg')
  print(res$gene[res$log2FoldChange < 0])
  
  #next
  ## Boxplot
  if (boxplot_TF){
    for (gene in res$gene){
      
      log2FoldChange = res$log2FoldChange[res$gene == gene]
      log2FoldChange = format(round(log2FoldChange, 2), nsmall = 2)
      filename <- paste0(path,paste0('HeatMap','.png'))
      
      folder = paste0(path,'BoxPlots/')
      dir.create(folder,recursive = TRUE)
      data = as.data.frame( DE_input@assays[["RNA"]]@data[gene,])
      colnames(data) = 'count'
      data$DE_ident = DE_input$DE_ident
      pathName = as.character(paste0(folder,'Boxplot_','Log2FC_',log2FoldChange,'_',gene,'.png') )
      png(file=pathName,width=900, height=600)
      title = paste0(gene,' Log2FC: ',log2FoldChange)
      plot = ggplot(data, aes(x=DE_ident, y=count)) + 
        geom_boxplot()
      plot = plot + ggtitle(title) +
        theme_classic()  + theme(
          title =element_text(size=24, color="black"),
          axis.title.x = element_text(color="black", size=24 ),
          axis.title.y = element_text(color="black", size=24),
          axis.text= element_text(color="black", size=24),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18),
          legend.position = 'none'
        )
      print(plot)
      dev.off()
    }
    next
  }

  
  ## Heatmap
  filename <- paste0(path,paste0('HeatMap_pos','.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res_pos),res = 100)
  
  plot = DoHeatmap(DE_input, features = res_pos$gene, group.by = 'DE_ident',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
  filename <- paste0(path,paste0('HeatMap_neg','.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res_neg),res = 100)
  
  plot = DoHeatmap(DE_input, features = res_neg$gene, group.by = 'DE_ident',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
  DE_input_ident1 = DE_input[,DE_input$DE_ident == ident1]
  DE_input_ident2 = DE_input[,DE_input$DE_ident == ident2]
  
  
  filename <- paste0(path,paste0('HeatMap','_patient_',ident1,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident1, features = res$gene, group.by = 'Patient Number',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = font_size))
  print(plot)
  dev.off()
  
  filename <- paste0(path,paste0('HeatMap','_patient_',ident2,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident2, features = res$gene, group.by = 'Patient Number',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  ###
  
  filename <- paste0(path,paste0('HeatMap','_Sample_',ident1,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident1, features = res$gene, group.by = 'Sample',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
  filename <- paste0(path,paste0('HeatMap','_Sample_',ident2,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident2, features = res$gene, group.by = 'Sample',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
  ###
  
  filename <- paste0(path,paste0('HeatMap','_dexa_',ident1,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident1, features = res$gene, group.by = 'Dexa or not',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  
  print(plot)
  dev.off()
  
  filename <- paste0(path,paste0('HeatMap','_dexa_',ident2,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident2, features = res$gene, group.by = 'Dexa or not',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  
  print(plot)
  dev.off()
  
  ###
  filename <- paste0(path,paste0('HeatMap','_kit','.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input, features = res$gene, group.by = 'kit',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  ## Plot sex
  
  
  filename <- paste0(path,paste0('HeatMap','_gender_',ident1,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident1, features = res$gene, group.by = 'Gender',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  
  print(plot)
  dev.off()
  
  filename <- paste0(path,paste0('HeatMap','_gender_',ident2,'.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input_ident2, features = res$gene, group.by = 'Gender',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  
  print(plot)
  dev.off()
  
  ###
 
  
  ##
  next
  # All time points
  ident1 = paste0('NBM ',celltype)
  ident2 = paste0('baseline ',celltype)
  ident3 = paste0('C9D1 ',celltype)
  ident4 = paste0('EOT ',celltype)
  
  DE_input = data_merge_run_label
  DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2,ident3,ident4)]
  
  DE_input$DE_ident = factor(DE_input$DE_ident,levels = c(ident1,ident2,ident3,ident4))
  
  filename <- paste0(path,paste0('HeatMap','_DE_ident_allTimes','.png'))
  png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input, features = res$gene, group.by = 'DE_ident',size = label_size)
  plot = plot + theme(
    axis.title.x = element_text(color="black", size=font_size ),
    axis.title.y = element_text(color="black", size=font_size),
    axis.text= element_text(color="black", size=font_size),
    legend.text=element_text(size=font_size),
    legend.title=element_text(size=font_size),
    text = element_text(size = 20))
  print(plot)
  dev.off()
  
  
}


#########################################
## Make a list of all DE genes and which timepoints they are DE
## Find genes that continue going up/down 
#########################################

colnames = c('padj_old','CellType','Time',"X","baseMean","log2FoldChange","lfcSE","stat", "pvalue" , "padj","name","gene")
allDE <- data.frame(matrix(ncol = length(colnames), nrow = 0))
colnames(allDE) <- colnames
celltype_list = sort(unique(as.character(Idents(data_merge_run_label))))


DE_type = 'DESeq2'
DE_type = 'EdgeR'

base = paste0(filepath_cluster,'DE/',DE_type,'/Patient/')



for (celltype in celltype_list){
  print(celltype)
  ident1 = paste0('baseline ',celltype)
  ident2 = paste0('C9D1 ',celltype)
  
  subfolder = paste0(ident1,' Vs ',ident2)
  
  path = paste0(base,'/baseline Vs C9D1/', subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',ident1,' Vs ',ident2,'.csv')
  
  if(file.exists(filename)){
    res = read.csv(filename)
    print('Exists baseline Vs C9D1')
    
    res$gene = res$X
    res$padj_old =  res$padj
    #res$padj  = res$FDR
    #res$padj = p.adjust(res$pvalue, method = 'BH', n = length(res$pvalue))
    
    names(res)[names(res) == "logFC"] <- "log2FoldChange"
    names(res)[names(res) == "FDR"] <- "padj"
    names(res)[names(res) == "PValue"] <- "pvalue"
    
    res = res[order(res$log2FoldChange),]
    res = res[res$padj < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl("^RP[SL]", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'baselineVsC9D1'
      allDE = rbind(allDE,res)
    }
    
    
  }
  
  # Get C9D1 Vs EOT
  ident1 = paste0('C9D1 ',celltype)
  ident2 = paste0('EOT ',celltype)
  
  subfolder = paste0( '/', ident1,' Vs ',ident2)
  
  path = paste0(base,'/C9D1 Vs EOT/', subfolder,'/')

  filename = paste0(path,'/DE_',DE_type,' ',ident1,' Vs ',ident2,'.csv')
  
  if(file.exists(filename)){
    res = read.csv(filename)
    print('Exists C9D1 Vs EOT')
    
    res$gene = res$X
    res$padj_old =  res$padj
    #res$padj = p.adjust(res$pvalue, method = 'BH', n = length(res$pvalue))
    
    names(res)[names(res) == "logFC"] <- "log2FoldChange"
    names(res)[names(res) == "FDR"] <- "padj"
    names(res)[names(res) == "PValue"] <- "pvalue"
    
    
    res = res[order(res$log2FoldChange),]
    res = res[res$padj < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl("^RP[SL]", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'C9D1VsEOT'
      allDE = rbind(allDE,res)
    }
    
    
  }
  
  # Get baseline Vs NBM
  ident1 = paste0('baseline ',celltype)
  ident2 = paste0('NBM ',celltype)
  
  subfolder = paste0( '/', ident1,' Vs ',ident2)
  
  path = paste0(base,'/baseline Vs NBM/', subfolder,'/')
  
  filename = paste0(path,'/DE_',DE_type,' ',ident1,' Vs ',ident2,'.csv')
  
  if(file.exists(filename)){
    res = read.csv(filename)
    print('Exists baseline Vs NBM')
    
    res$gene = res$X
    res$padj_old =  res$padj
    #res$padj = p.adjust(res$pvalue, method = 'BH', n = length(res$pvalue))
    
    names(res)[names(res) == "logFC"] <- "log2FoldChange"
    names(res)[names(res) == "FDR"] <- "padj"
    names(res)[names(res) == "PValue"] <- "pvalue"
    
    
    res = res[order(res$log2FoldChange),]
    res = res[res$padj < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl("^RP[SL]", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'baselineVsNBM'
      allDE = rbind(allDE,res)
    }
    
    
  }
  
  # Get baseline GR Vs baseline PR
  ident1 = paste0('baseline GR ',celltype)
  ident2 = paste0('baseline PR ',celltype)
  
  subfolder = paste0( '/', ident1,' Vs ',ident2)
  
  path = paste0(base,'/baseline GR Vs baseline PR/', subfolder,'/')
  
  filename = paste0(path,'/DE_',DE_type,' ',ident1,' Vs ',ident2,'.csv')
  
  if(file.exists(filename)){
    res = read.csv(filename)
    print('Exists baseline GR Vs baseline PR')
    
    res$gene = res$X
    res$padj_old =  res$padj
    #res$padj = p.adjust(res$pvalue, method = 'BH', n = length(res$pvalue))
    
    names(res)[names(res) == "logFC"] <- "log2FoldChange"
    names(res)[names(res) == "FDR"] <- "padj"
    names(res)[names(res) == "PValue"] <- "pvalue"
    
    
    res = res[order(res$log2FoldChange),]
    res = res[res$padj < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl("^RP[SL]", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'baselineGRVsbaselinePR'
      allDE = rbind(allDE,res)
    }
    
    
  }
}


############################################################################
## Get genes that are both different in baseline Vs C9D1 and C9D1 vs EOT
#############################################################################

celltype_list = unique(allDE$CellType)

for (celltype in celltype_list){
  sig_pos = c()
  sig_neg = c()
  allDE_celltype = allDE[allDE$CellType == celltype,]
  
  gene_list = unique(allDE_celltype$gene)
  
  
  allDE_pos = allDE_celltype[allDE_celltype$log2FoldChange > 0,]
  allDE_neg = allDE_celltype[allDE_celltype$log2FoldChange < 0,]
  for (gene in gene_list){
    
    time_list_pos = unique(allDE_pos$time[allDE_pos$gene == gene])
    time_list_neg = unique(allDE_neg$time[allDE_neg$gene == gene])
    
    if( c('baselineVsC9D1','C9D1VsEOT') %in% time_list_pos  ){
      sig_pos = c(sig_pos,gene)
    }
    
    if( c('baselineVsC9D1','C9D1VsEOT') %in% time_list_neg){
      sig_neg = c(sig_neg,gene)
    }
    
    
  }
  
  if (length(sig_pos) > 0 | length(sig_neg) > 0){
    ident1 = paste0('NBM ',celltype)
    ident2 = paste0('baseline ',celltype)
    ident3 = paste0('C9D1 ',celltype)
    ident4 = paste0('EOT ',celltype)
    
    DE_input = data_merge_run_label
    DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
    DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2,ident3,ident4)]
    
    DE_input$DE_ident = factor(DE_input$DE_ident,levels = c(ident1,ident2,ident3,ident4))
  }
  
  if (length(sig_pos) > 0){

    print(celltype)
    print('sig_pos')
    print(sig_pos)
    print(' ')
    
    filename <- paste0(base,celltype,'_ContinuousPosHeatMap','.png')
    png(file=filename,width=2000, height=gene_height*length(sig_pos),res = 100)
    
    plot = DoHeatmap(DE_input, features = sig_pos, group.by = 'DE_ident',size = label_size)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=font_size ),
      axis.title.y = element_text(color="black", size=font_size),
      axis.text= element_text(color="black", size=font_size),
      legend.text=element_text(size=font_size),
      legend.title=element_text(size=font_size),
      text = element_text(size = 20))
    print(plot)
    dev.off()
  }
  
  if (length(sig_neg) > 0){
    print(celltype)
    print('sig_neg')
    print(sig_neg)
    print(' ')
    
    
    
    filename <- paste0(base,celltype,'_ContinuousNegHeatMap','.png')
    png(file=filename,width=2000, height=gene_height*length(sig_neg),res = 100)
    
    plot = DoHeatmap(DE_input, features = sig_neg, group.by = 'DE_ident',size = label_size)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=font_size ),
      axis.title.y = element_text(color="black", size=font_size),
      axis.text= element_text(color="black", size=font_size),
      legend.text=element_text(size=font_size),
      legend.title=element_text(size=font_size),
      text = element_text(size = 20))
    print(plot)
    dev.off()
  }
  
 

}

############################################################################
## Get genes that are both different in baseline Vs NBM and baseline vs C9D1
#############################################################################

for (celltype in celltype_list){
  sig_pos = c()
  sig_neg = c()
  
  allDE_celltype = allDE[allDE$CellType == celltype,]
  
  gene_list = unique(allDE_celltype$gene)
  
  for (gene in gene_list){
    
    time_list = unique(allDE_celltype$time[allDE_celltype$gene == gene])
    #print(time_list)
    # Get FC between baseline and C9
    log2FoldChange = allDE_celltype$log2FoldChange[allDE_celltype$gene == gene & 
                                                     allDE_celltype$time == 'baselineVsC9D1']
    if (length(log2FoldChange) == 0){
      log2FoldChange = 0
    }
    
    if( c('baselineVsNBM','baselineVsC9D1') %in% time_list  &  log2FoldChange > 0 ){
      sig_pos = c(sig_pos,gene)
    }
    
    if( c('baselineVsNBM','baselineVsC9D1') %in% time_list  &  log2FoldChange < 0 ){
      sig_neg = c(sig_neg,gene)
    }
    
    
  }
  
  ident1 = paste0('NBM ',celltype)
  ident2 = paste0('baseline ',celltype)
  ident3 = paste0('C9D1 ',celltype)
  ident4 = paste0('EOT ',celltype)
  
  DE_input = data_merge_run_label
  DE_input$DE_ident = paste0(DE_input$Treatment, ' ', Idents(DE_input))
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2,ident3,ident4)]
  
  DE_input$DE_ident = factor(DE_input$DE_ident,levels = c(ident1,ident2,ident3,ident4))
  
  if (length(sig_pos) > 0){
    
    print(celltype)
    print('sig_pos')
    print(sig_pos)
    print(' ')
    
    filename <- paste0(base,celltype,'_baselineVsNBMAndBaselineVsC9D1_pos','.png')
    png(file=filename,width=2000, height=gene_height*length(sig_pos) + 400,res = 100)
    
    plot = DoHeatmap(DE_input, features = sig_pos, group.by = 'DE_ident',size = label_size)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=font_size ),
      axis.title.y = element_text(color="black", size=font_size),
      axis.text= element_text(color="black", size=font_size),
      legend.text=element_text(size=font_size),
      legend.title=element_text(size=font_size),
      text = element_text(size = 20))
    print(plot)
    dev.off()
  }
  
  if (length(sig_neg) > 0){
    
    print(celltype)
    print('sig_neg')
    print(sig_neg)
    print(' ')
    
    filename <- paste0(base,celltype,'_baselineVsNBMAndBaselineVsC9D1_neg','.png')
    png(file=filename,width=2000, height=gene_height*length(sig_neg) + 400,res = 100)
    
    plot = DoHeatmap(DE_input, features = sig_neg, group.by = 'DE_ident',size = label_size)
    plot = plot + theme(
      axis.title.x = element_text(color="black", size=font_size ),
      axis.title.y = element_text(color="black", size=font_size),
      axis.text= element_text(color="black", size=font_size),
      legend.text=element_text(size=font_size),
      legend.title=element_text(size=font_size),
      text = element_text(size = 20))
    print(plot)
    dev.off()
  }
  
}

##
## Look at individual genes
##
celltype = 'Naive CD8+ T-cell'
gene = 'CD48'
allDE_celltype = allDE[allDE$CellType == celltype,]
tmp = allDE_celltype[allDE_celltype$gene == gene & allDE_celltype$time == 'baselineVsC9D1',]
tmp = allDE_celltype[allDE_celltype$gene == gene ,]

tmp
#######################################
## Analyze Patient data
## For each cell type comparision, load all the available patients
## find gene that are up/down in 3 or more patients
#####################################3

base = paste0(filepath_cluster,'DE/DESeq2/Individual Patient/')


patient_list = sort(unique(data_merge_run_label$`Patient Number`))

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'dMono','Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype_list=  c('eTreg')
celltype = 'CD4+ TCM'
patient = 1

colnames = c('Patient','CellType',"X","baseMean","log2FoldChange","lfcSE","stat", "pvalue" , "padj","name","gene")
allDE <- data.frame(matrix(ncol = length(colnames), nrow = 0))
colnames(allDE) <- colnames


for (patient in patient_list){
  for (celltype in celltype_list){
    print(celltype)
    ident1 = paste0('baseline ',celltype)
    ident2 = paste0('C9D1 ',celltype)
    
    subfolder = paste0('Patient',patient, '/', ident1,' Vs ',ident2)
    
    path = paste0(base, subfolder,'/')
    filename = paste0(path,'/DE_DESeq2 ',ident1,' Vs ',ident2,'.csv')
    
    if(file.exists(filename)){
      res = read.csv(filename)
      print('Exists')
    
      res$gene = res$X
      res$padj = p.adjust(res$pvalue, method = 'BH', n = length(res$pvalue))
      
      
      res = res[order(res$log2FoldChange),]
      res = res[res$padj < 0.05,]
      
      res = res[!grepl("MT-", res$gene),]
      res = res[!grepl("^RP[SL]", res$gene),]
      
      res = res[rowSums(is.na(res)) != ncol(res),]
      
      if (nrow(res) > 0){
        res$Patient = patient
        res$CellType = celltype
        allDE = rbind(allDE,res)
      }
      
      
    }
  }
}




for (celltype in celltype_list){
  print(celltype)
  sig_pos = c()
  sig_neg = c()
  allDE_celltype = allDE[allDE$CellType == celltype,]

  gene_list = unique(allDE_celltype$gene)
  
  
  allDE_pos = allDE_celltype[allDE_celltype$log2FoldChange > 0,]
  allDE_neg = allDE_celltype[allDE_celltype$log2FoldChange < 0,]
  
  for (gene in gene_list){
    
    patient_list_pos = unique(allDE_pos$Patient[allDE_pos$gene == gene])
    patient_list_neg = unique(allDE_neg$Patient[allDE_neg$gene == gene])
    
    if( length(patient_list_pos ) > 2){
      sig_pos = c(sig_pos,gene)
    }
    
    if( length(patient_list_neg) > 2){
      sig_neg = c(sig_neg,gene)
    }
   

  }
  
  if (length(sig_pos) > 0){
    print('sig_pos')
    print(sig_pos)
    print(' ')
  }
  
  if (length(sig_neg) > 0){
    print('sig_neg')
    print(sig_neg)
    print(' ')
  }
}


#######################
## Check which genes overlap between time and PR Vs GR
#######################

celltype_list=  c('Th2')
celltype = celltype_list[1]
DE_type = 'EdgeR'
time1 = 'baseline'
time2 = 'C9D1'

time = 'EOT'

for (celltype in celltype_list){
  print(celltype)
  ident1 = paste0(time1,' ',celltype)
  ident2 = paste0(time2,' ',celltype)
  
  base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/', time1,' Vs ',time2,'/')
  
  subfolder = paste0(ident1,' Vs ',ident2)

  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'.csv')
  
  if(!file.exists(filename)){
    next
  }
  
  res_time = read.csv(filename)
  res_time$gene = res_time$X
  names(res_time)[names(res_time) == "logFC"] <- "log2FoldChange"
  names(res_time)[names(res_time) == "FDR"] <- "padj"
  names(res_time)[names(res_time) == "PValue"] <- "pvalue"
  
  res_time = res_time[order(res_time$log2FoldChange),]
  res_time = res_time[res_time$padj < 0.05,]
  
  res_time = res_time[!grepl("MT-", res_time$gene),]
  res_time = res_time[!grepl("^RP[SL]", res_time$gene),]
  
  res_time = res_time[rowSums(is.na(res_time)) != ncol(res_time),]
  
  
  ident1 = paste0(time,' GR ',celltype)
  ident2 = paste0(time,' PR ',celltype)
  
  subfolder = paste0(ident1,' Vs ',ident2)
  
  base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/', time,' GR',' Vs ',time,' PR','/')
  
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'.csv')
  
  if(!file.exists(filename)){
    next
  }
  
  res_response = read.csv(filename)
  res_response$gene = res_response$X
  names(res_response)[names(res_response) == "logFC"] <- "log2FoldChange"
  names(res_response)[names(res_response) == "FDR"] <- "padj"
  names(res_response)[names(res_response) == "PValue"] <- "pvalue"
  
  res_response = res_response[order(res_response$log2FoldChange),]
  res_response = res_response[res_response$padj < 0.05,]
  
  res_response = res_response[!grepl("MT-", res_response$gene),]
  res_response = res_response[!grepl("^RP[SL]", res_response$gene),]
  
  res_response = res_response[rowSums(is.na(res_response)) != ncol(res_response),]
  
  max_gene = 50
  res_time_short = res_time
  res_time_short$gene = as.character(res_time_short$gene)
  res_pos = res_time_short[res_time_short$log2FoldChange > 0,]
  res_neg = res_time_short[res_time_short$log2FoldChange < 0,]
  if (nrow(res_pos) > max_gene){
    
    res_pos = res_pos[1:max_gene,]
    res_pos = res_pos[rowSums(is.na(res_pos)) != ncol(res_pos),]
    
  }
  if (nrow(res_neg) > max_gene){
    
    res_neg = res_neg[( nrow(res_time_short) - max_gene):nrow(res_time_short),]
    res_neg = res_neg[rowSums(is.na(res_pos)) != ncol(res_pos),]
    
  }
  res_time_short = res_time_short[res_time_short$gene %in% c(res_pos$gene, res_neg$gene),]
  res_time_short = unique(res_time_short)
  res_time_short = res_time_short[rowSums(is.na(res_time_short)) != ncol(res_time_short),]
  
  res_time = res_time_short
  print('Positive intersection')
  print(res_time_short$gene[res_time_short$log2FoldChange > 0],)
  print(intersect(res_time$gene[res_time$log2FoldChange > 0],res_response$gene))
  gene_list = intersect(res_time$gene[res_time$log2FoldChange > 0],res_response$gene)
  print(res_response[res_response$gene %in% gene_list,])
  
  print('Negative intersection')
  print(res_time_short$gene[res_time_short$log2FoldChange < 0],)
  print(intersect(res_time$gene[res_time$log2FoldChange < 0],res_response$gene))
  gene_list = intersect(res_time$gene[res_time$log2FoldChange < 0],res_response$gene)
  print(res_response[res_response$gene %in% gene_list,])
}

