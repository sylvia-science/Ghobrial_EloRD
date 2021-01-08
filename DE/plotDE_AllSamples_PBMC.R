

data_harmony_run_label_subset = ScaleData(data_harmony_run_label, features = rownames(data_harmony_run_label))
#data_harmony_run_label = renameCells(data_harmony_run_label,idents = c('cDC1','cDC2'),newident = 'DC')
#data_harmony_run_label = renameCells(data_harmony_run_label,idents = c('CD14+CD16+ Mono'),newident = 'CD16+ Mono')

remove_list = c('HSC','Pro Erythrocyte','Remove','Plasma Cell','GMPC','Pro B Cell','Pre B Cell','CMPC','MDPC')
data_harmony_run_label_subset = data_harmony_run_label_subset[,!(Idents(data_harmony_run_label_subset) %in% remove_list)]

DE_input = data_harmony_run_label_subset[,Idents(data_harmony_run_label_subset) == 'T-cell']
DE_input = DE_input[,DE_input$Treatment %in% c('baseline','C9D1')]
DE_input$DE_ident = paste0(DE_input$Treatment,' ', DE_input$`Sample Type` , ' ', Idents(DE_input))
table(DE_input$`Patient Number`, DE_input$DE_ident )

patientByTreamtment = table(DE_input$`Patient Number`, paste0(DE_input$Treatment,' ', DE_input$`Sample Type`))

write.csv(patientByTreamtment, file = paste0(base,'patientByTreamtment.csv'),row.names=TRUE)

library(fgsea)
pathways.hallmark <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/h.all.v7.2.symbols.gmt')
pathways.c2 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c2.all.v7.2.symbols.gmt')
pathways.c5 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c5.all.v7.2.symbols.gmt')
pathways.c7 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c7.all.v7.2.symbols.gmt')
pathways.c8 <- gmtPathways('/home/sujwary/Desktop/scRNA/Data/GSEA/c8.all.v7.2.symbols.gmt')
pathways = c(pathways.hallmark,pathways.c2,pathways.c5, pathways.c7,pathways.c8)


font_size= 32

label_size = 12

gene_height = 70

celltype_list = sort(unique(as.character(Idents(data_harmony_run_label_subset))))
for (celltype in celltype_list){
  print(celltype)
}


DE_type = 'DESeq2'
DE_type = 'EdgeR'

time1 = 'baseline BM'
time2 = 'baseline PB'

base = paste0(filepath_cluster,'/DE/',DE_type,'/')
#base = paste0(filepath_cluster,'/DE/',DE_type,'/Patient/C9D1 Dexa Vs C9D1 No Dexa/')



celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','Plasma Cell')


celltype_list=  c('Intermediate CD4+ T-cell')
celltype = celltype_list[1]

celltype = unique()

boxplot_TF = F

cell_features = getCellMarkers('/home/sujwary/Desktop/scRNA/')

cell_features$Cell[cell_features$Cell == 'tcell_general'] = 'T-cell'
cell_features$Cell[cell_features$Cell == 'tcell_cytotoxic'] = 'CD8+ T-cell'
cell_features$Cell[cell_features$Cell == 'monocyte_CD14'] = 'CD14+ Mono'
cell_features$Cell[cell_features$Cell == 'nk_cell'] = 'NK'
cell_features$Cell[cell_features$Cell == 'monocyte_FCGR3A'] = 'CD16+ Mono'
cell_features$Cell[cell_features$Cell == 'bcell_immature'] = 'B-cell'
cell_features$Cell[cell_features$Cell == 'bcell_mature'] = 'B-cell'


celltype = celltype_list[1]

features = read.csv('/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Cluster/PCA40/res3/Features_label.csv')

for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0(time1,' ',celltype)
  ident2 = paste0(time2,' ',celltype)
  

  
  subfolder = paste0(ident1,' Vs ',ident2)
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'.csv')
  print(file.exists(filename))
  if(!file.exists(filename)){
    next
  }
  
  
  DE_input = data_harmony_run_label_subset
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Sample_type = ''
  DE_input$Sample_type[DE_input$'Sample Type' == "PBMC"] = 'PB'
  DE_input$Sample_type[DE_input$'Sample Type' == "Bone Marrow"] = 'BM'
  
  DE_input$DE_ident = paste0(DE_input$Treatment,' ', DE_input$Sample_type, ' ', Idents(DE_input))
  
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
  
  
  res_orig = read.csv(filename)
  res_orig$gene = res_orig$X
  names(res_orig)[names(res_orig) == "logFC"] <- "log2FoldChange"
  names(res_orig)[names(res_orig) == "FDR"] <- "padj"
  names(res_orig)[names(res_orig) == "PValue"] <- "pvalue"
  
  res_orig = res_orig[!grepl("MT-", res_orig$gene),]
  res_orig = res_orig[!grepl("^RP[SL]", res_orig$gene),]
  
  res = res_orig[order(res_orig$log2FoldChange),]
  res = res[res$padj < 0.05,]
  
  
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/DE_',DE_type,' ',subfolder,'_clean','.csv')
  
  write.csv(res, file = filename,row.names=TRUE)
  
  #next
  

  
  res = res[rowSums(is.na(res)) != ncol(res),]
  

  
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
  
  res_pos = res_pos[order(-res_pos$log2FoldChange),]
  
  
  if (nrow(res_pos) > max_gene){
    
    res_pos = res_pos[1:max_gene,]
    res_pos = res_pos[rowSums(is.na(res_pos)) != ncol(res_pos),]
    
  }
  if (nrow(res_neg) > max_gene){
    
    res_neg = res_neg[1:max_gene,]
    res_neg = res_neg[rowSums(is.na(res_neg)) != ncol(res_neg),]
    
  }
  
  
  #print('pos')
  #print(res$gene[res$log2FoldChange > 0])
  
  #print('neg')
  #print(res$gene[res$log2FoldChange < 0])
  
  folder = paste0(filepath_cluster,'/DE/',DE_type,'/',subfolder,'/')
  
  res_pos$gene %in% rownames(DE_input@assays[["RNA"]]@scale.data)
  if(nrow(res_pos) > 0){
    filename <- paste0(folder,paste0('HeatMap_pos','.png'))
    png(file=filename,width=2000, height=gene_height*nrow(res_pos),res = 100)
    
    plot = DoHeatmap(DE_input,slot = 'scale.data', features = res_pos$gene, group.by = 'DE_ident',size = label_size)
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
  
  if(nrow(res_neg) > 0){
    
    filename <- paste0(folder,paste0('HeatMap_neg','.png'))
    png(file=filename,width=2000, height=gene_height*nrow(res_neg),res = 100)
    
    plot = DoHeatmap(DE_input,slot = 'scale.data', features = res_neg$gene, group.by = 'DE_ident',size = label_size)
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
  
  gene_list = cell_features$Markers[cell_features$Cell == celltype]
  if (length(gene_list) > 1){
    gene_list = paste(gene_list,collapse="")
  }
  
  gene_list =   unlist(strsplit(gene_list, ",")) 

  gene_list = trimws(gene_list, which = c("both"), whitespace = " \t\n\r\v\f")
  gene_list = gsub(" ", "", gene_list)

  gene_list = unique(gene_list)

  res_marker = res_orig[res_orig$gene %in% gene_list,]

  colnames(res_marker) <- make.unique(names(res_marker))
  library("grid")
  library(ggrepel)
  path = paste0(base, subfolder,'/')
  filename = paste0(path,'/cellMarkers',celltype,'.csv')
  write.csv(gene_list, file = filename,row.names=F)
  
  
  if (nrow(res_marker) > 0){
    print(res_marker)
    ymax = ceiling(max(-log10(res_marker$padj) )) + 1
    ymin= floor(min(-log10(res_marker$padj) )) -1
    xmax = 3 #ceiling(max(res_marker$log2FoldChange) ) + 1
    xmin = -3 #floor(min(res_marker$log2FoldChange) ) - 1

    crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
    g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
    
    pdf(paste0(folder,'Composition Volcano Plot_',ident1,'Vs',ident2,'.pdf'))
    
    plot = ggplot(res_marker,aes(x=log2FoldChange, y=-log10(padj)), size=4) +
      annotation_custom(g, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
      geom_point() +
      xlim(xmin,xmax) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")
    plot = plot + xlab("Log2 fold-change")+ylab("-log10 p-value") + 
      theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    plot = plot + geom_text_repel(aes(x=log2FoldChange, y=-log10(padj),label=gene))+
      geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
    
    print(plot)
    dev.off()
  
  }
  
  # Plot top DE features of each cell type
  
  gene_list = as.character(features$gene[features$cluster == celltype])
  gene_list = gene_list[!grepl("MT-", gene_list)]
  gene_list = gene_list[!grepl("^RP[SL]", gene_list)]
  
  gene_list = gene_list[1:20]

  res_marker = res_orig[res_orig$gene %in% gene_list,]
  colnames(res_marker) <- make.unique(names(res_marker))
  if (nrow(res_marker) > 0){
    print(res_marker)
    ymax = ceiling(max(-log10(res_marker$padj) )) + 1
    ymin= floor(min(-log10(res_marker$padj) )) -1
    xmax = 3 #ceiling(max(res_marker$log2FoldChange) ) + 1
    xmin = -3 #floor(min(res_marker$log2FoldChange) ) - 1
    
    crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
    g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
    
    pdf(paste0(folder,'Volcano Plot_DEFeatures_',ident1,'Vs',ident2,'.pdf'))
    
    plot = ggplot(res_marker,aes(x=log2FoldChange, y=-log10(padj)), size=4) +
      annotation_custom(g, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
      geom_point() +
      xlim(xmin,xmax) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")
    plot = plot + xlab("Log2 fold-change")+ylab("-log10 p-value") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    plot = plot + geom_text_repel(aes(x=log2FoldChange, y=-log10(padj),label=gene))+
      geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
    
    print(plot)
    dev.off()
    
  }
  
  #####
  
  #next
  ## Boxplot
  if (boxplot_TF){
    for (gene in res$gene){
      
      log2FoldChange = res$log2FoldChange[res$gene == gene]
      log2FoldChange = format(round(log2FoldChange, 2), nsmall = 2)

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
  
  
  # ## Heatmap for overlapping marker genes
  # filename <- paste0(path,paste0('HeatMap','.png'))
  # png(file=filename,width=2000, height=gene_height*nrow(res),res = 100)
  # 
  # plot = DoHeatmap(DE_input, features = res$gene, group.by = 'DE_ident',size = label_size)
  # plot = plot + theme(
  #   axis.title.x = element_text(color="black", size=font_size ),
  #   axis.title.y = element_text(color="black", size=font_size),
  #   axis.text= element_text(color="black", size=font_size),
  #   legend.text=element_text(size=font_size),
  #   legend.title=element_text(size=font_size),
  #   text = element_text(size = 20))
  # print(plot)
  # dev.off()
  

  
  next
  
    
  
}



patientByGroup =  table(data_harmony_run_label_subset$`Patient Number`,
                        paste0(data_harmony_run_label_subset$Treatment,' ', data_harmony_run_label_subset$'Sample Type'))

patientByGroup = as.data.frame(patientByGroup > 0)
patient_list = rownames(patientByGroup)[patientByGroup$`baseline Bone Marrow` & patientByGroup$`baseline PBMC`]

#patient_list = sort(unique(data_harmony_run_label$`Patient Number`))
#patient_list

patient = patient_list[1]


for (patient in patient_list){
  print(patient)
  for (celltype in celltype_list){
    print(celltype)
    #celltype = 'eTreg'
    ident1 = paste0('baseline BM ',celltype)
    ident2 = paste0('baseline PB ',celltype)
    
    subfolder = paste0('Patient',patient,'/',ident1,' Vs ',ident2)
    
    
    folder = paste0(filepath_cluster,'/DE/EdgeR/',subfolder,'/')
    
    filename = paste0(folder,'Patient',patient,'_clean','.csv')
    if(!file.exists(filename)){
      next
    }
    
    #ident1 = paste0('CD16+ Mono')
    #ident2 = paste0('CD14+ Mono')
    
    
    DE_input = data_harmony_run_label_subset
    print( unique(DE_input$`Sample Type`))
    DE_input = DE_input[,DE_input$`Patient Number` == patient]
    print( unique(DE_input$`Sample Type`))
    DE_input$Sample_type = ''
    DE_input$Sample_type[DE_input$'Sample Type' == "PBMC"] = 'PB'
    DE_input$Sample_type[DE_input$'Sample Type' == "Bone Marrow"] = 'BM'
    
    
    table(DE_input$Treatment, DE_input$Sample_type )
    DE_input$DE_ident = paste0(DE_input$Treatment,' ', DE_input$Sample_type, ' ', Idents(DE_input))
    
    print(    table(DE_input$`Patient Number`,paste0(DE_input$Treatment,' ', DE_input$'Sample Type')))
    
    #DE_input = DE_input[,Idents(DE_input) == celltype]
    ncol(DE_input)
    unique(DE_input$DE_ident)
    
    
    if (!( c(ident1,ident2) %in% DE_input$DE_ident)){
      print('Idents do not exist')
      next
    }
    DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
    ncol(DE_input)
    unique(DE_input$DE_ident)
    
    if (ncol(DE_input) < 100){
      print('Not enough cells')
      next
    }
   
   

    res_orig = read.csv(filename)
    res_orig$gene = res_orig$X
    names(res_orig)[names(res_orig) == "logFC"] <- "log2FoldChange"
    names(res_orig)[names(res_orig) == "FDR"] <- "padj"
    names(res_orig)[names(res_orig) == "PValue"] <- "pvalue"
    
    res_orig = res_orig[!grepl("MT-", res_orig$gene),]
    res_orig = res_orig[!grepl("^RP[SL]", res_orig$gene),]
    
    res = res_orig[order(res_orig$log2FoldChange),]
    
    
    res$pvalue = res$PValue
    res$log2FoldChange = res$logFC

    res = res[order(res$log2FoldChange),]
    res = res[res$pvalue < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl(" ?RP\\w+ ?", res$gene),]
    
    max_gene = 50
    res$gene = as.character(res$gene)
    res_pos = res[res$log2FoldChange > 0,]
    res_neg = res[res$log2FoldChange < 0,]
    
    res_pos = res_pos[order(-res_pos$log2FoldChange),]
    
    
    if (nrow(res_pos) > max_gene){
      
      res_pos = res_pos[1:max_gene,]
      res_pos = res_pos[rowSums(is.na(res_pos)) != ncol(res_pos),]
      
    }
    if (nrow(res_neg) > max_gene){
      
      res_neg = res_neg[1:max_gene,]
      res_neg = res_neg[rowSums(is.na(res_neg)) != ncol(res_neg),]
      
    }
    
    
    #print('pos')
    #print(res$gene[res$log2FoldChange > 0])
    
    #print('neg')
    #print(res$gene[res$log2FoldChange < 0])
    
    folder = paste0(filepath_cluster,'/DE/',DE_type,'/',subfolder,'/')
    
    res_pos$gene %in% rownames(DE_input@assays[["RNA"]]@scale.data)
    
    if(nrow(res_pos) > 0){
      filename <- paste0(folder,paste0('HeatMap_pos','.png'))
      png(file=filename,width=2000, height=gene_height*nrow(res_pos),res = 100)
      
      plot = DoHeatmap(DE_input,slot = 'scale.data', features = res_pos$gene, group.by = 'DE_ident',size = label_size)
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
    
    if(nrow(res_neg) > 0){
      
      filename <- paste0(folder,paste0('HeatMap_neg','.png'))
      png(file=filename,width=2000, height=gene_height*nrow(res_neg),res = 100)
      
      plot = DoHeatmap(DE_input,slot = 'scale.data', features = res_neg$gene, group.by = 'DE_ident',size = label_size)
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
    ###################
    gene_list = cell_features$Markers[cell_features$Cell == celltype]
    if (length(gene_list) > 1){
      gene_list = paste(gene_list,collapse="")
    }
    
    gene_list =   unlist(strsplit(gene_list, ",")) 
    
    gene_list = trimws(gene_list, which = c("both"), whitespace = " \t\n\r\v\f")
    gene_list = gsub(" ", "", gene_list)
    
    gene_list = unique(gene_list)
    
    res_marker = res_orig[res_orig$gene %in% gene_list,]
    
    colnames(res_marker) <- make.unique(names(res_marker))
    library("grid")
    library(ggrepel)
    path = paste0(base, subfolder,'/')
    filename = paste0(folder,'/cellMarkers',celltype,'.csv')
    write.csv(gene_list, file = filename,row.names=F)
    
    
    if (nrow(res_marker) > 0){
      print(res_marker)
      ymax = ceiling(max(-log10(res_marker$padj) )) + 1
      ymin= floor(min(-log10(res_marker$padj) )) -1
      xmax = 3 #ceiling(max(res_marker$log2FoldChange) ) + 1
      xmin = -3 #floor(min(res_marker$log2FoldChange) ) - 1
      
      crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
      g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
      
      pdf(paste0(folder,'Composition Volcano Plot_',ident1,'Vs',ident2,'.pdf'))
      
      plot = ggplot(res_marker,aes(x=log2FoldChange, y=-log10(padj)), size=4) +
        annotation_custom(g, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
        geom_point() +
        xlim(xmin,xmax) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")
      plot = plot + xlab("Log2 fold-change")+ylab("-log10 p-value") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
      plot = plot + geom_text_repel(aes(x=log2FoldChange, y=-log10(padj),label=gene))+
        geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
      
      print(plot)
      dev.off()
      
    }
    
    # Plot top DE features of each cell type
    
    gene_list = as.character(features$gene[features$cluster == celltype])
    gene_list = gene_list[!grepl("MT-", gene_list)]
    gene_list = gene_list[!grepl("^RP[SL]", gene_list)]
    
    gene_list = gene_list[1:20]
    
    res_marker = res_orig[res_orig$gene %in% gene_list,]
    colnames(res_marker) <- make.unique(names(res_marker))
    if (nrow(res_marker) > 0){
      print(res_marker)
      ymax = ceiling(max(-log10(res_marker$padj) )) + 1
      ymin= floor(min(-log10(res_marker$padj) )) -1
      xmax = 3 #ceiling(max(res_marker$log2FoldChange) ) + 1
      xmin = -3 #floor(min(res_marker$log2FoldChange) ) - 1
      
      crange <- t(matrix(c("#009BF4","#EAEAEA","#FC5A5A","#EAEAEA","#EAEAEA","#EAEAEA"),ncol=2))
      g <- rasterGrob(crange, width=unit(1,"npc"), height = unit(1,"npc"),interpolate = TRUE)
      
      pdf(paste0(folder,'Volcano Plot_DEFeatures_',ident1,'Vs',ident2,'.pdf'))
      
      plot = ggplot(res_marker,aes(x=log2FoldChange, y=-log10(padj)), size=4) +
        annotation_custom(g, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax) +
        geom_point() +
        xlim(xmin,xmax) + geom_vline(xintercept = 0, linetype="dashed", alpha=0.5, color="black")
      plot = plot + xlab("Log2 fold-change")+ylab("-log10 p-value") + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
      plot = plot + geom_text_repel(aes(x=log2FoldChange, y=-log10(padj),label=gene))+
        geom_hline(yintercept = 1.3, linetype = "dashed", alpha = 0.5)
      
      print(plot)
      dev.off()
      
    }
    
    
   
  }
}

