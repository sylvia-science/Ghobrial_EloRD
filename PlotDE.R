
data_merge_run_label = ScaleData(data_merge_run_label, features = rownames(data_merge_run_label))

#data_merge_run_label <- NormalizeData(data_merge_run_label, normalization.method = "LogNormalize", scale.factor = 10000)

data_merge_run_label = renameCells(data_merge_run_label,idents = c('cDC1','cDC2'),newident = 'DC')
data_merge_run_label = renameCells(data_merge_run_label,idents = c('CD14+CD16+ Mono'),newident = 'CD16+ Mono')

library(fgsea)
pathways.hallmark <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/h.all.v7.2.symbols.gmt')
pathways.c2 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c2.all.v7.2.symbols.gmt')
pathways.c5 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c5.all.v7.2.symbols.gmt')
pathways.c7 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c7.all.v7.2.symbols.gmt')
pathways.c8 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c8.all.v7.2.symbols.gmt')
pathways = c(pathways.hallmark,pathways.c2,pathways.c5, pathways.c7,pathways.c8)


font_size= 32

label_size = 12

gene_height = 50

celltype_list = sort(unique(as.character(Idents(data_merge_run_label))))
for (celltype in celltype_list){
  print(celltype)
}


DE_type = 'DESeq2'
DE_type = 'EdgeR'

time1 = 'baseline'
time2 = 'C9D1'

base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/', time1,' Vs ',time2,'/')
#base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/C9D1 Dexa Vs C9D1 No Dexa/')



celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','Plasma Cell')


celltype_list=  c('Naive CD4+ T-cell')
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
  
  
  
  DE_input = data_merge_run_label
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
  
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl(" ?RP\\w+ ?", res$gene),]
  
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
  len = nrow(res)
  if (nrow(res) > max_gene*2){
    res = res[c(1:max_gene,(len - max_gene + 1):len),]
    res = res[rowSums(is.na(res)) != ncol(res),]
    res = unique(res)
  }
  res = res[rowSums(is.na(res)) != ncol(res),]
  res$gene = as.character(res$gene)
  
  print('pos')
  print(res$gene[res$log2FoldChange > 0])
  
  print('neg')
  print(res$gene[res$log2FoldChange < 0])
  
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
  filename <- paste0(path,paste0('HeatMap','.png'))
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
    res = res[!grepl(" ?RP\\w+ ?", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'baselineVsC9D1'
      allDE = rbind(allDE,res)
    }
    
    
  }
  
  
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
    res = res[!grepl(" ?RP\\w+ ?", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'C9D1VsEOT'
      allDE = rbind(allDE,res)
    }
    
    
  }
  
  
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
    res = res[!grepl(" ?RP\\w+ ?", res$gene),]
    
    res = res[rowSums(is.na(res)) != ncol(res),]
    
    if (nrow(res) > 0){
      res <- res[, !duplicated(colnames(res))]
      
      res$CellType = celltype
      res$time = 'baselineVsNBM'
      allDE = rbind(allDE,res)
    }
    
    
  }
}


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
      res = res[!grepl(" ?RP\\w+ ?", res$gene),]
      
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
