source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')
library(robustbase)
library(DESeq2)
library(EnrichmentBrowser)

font_size= 32

label_size = 12

gene_height = 70

data_harmony_run_label = ScaleData(data_harmony_run_label, features = rownames(data_harmony_run_label))
#data_harmony_run_label = renameCells(data_harmony_run_label,idents = c('cDC1','cDC2'),newident = 'DC')
#data_harmony_run_label = renameCells(data_harmony_run_label,idents = c('CD14+CD16+ Mono'),newident = 'CD16+ Mono')

remove_list = c('HSC','Pro Erythrocyte','Remove','Plasma-cell','GMPC','Pro B-cell','Pre B-cell','CMPC','MDPC',
                '13','30','41','49')
data_harmony_run_label_subset = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]


plot = DimPlot(data_harmony_run_label,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)
plot = DimPlot(data_harmony_run_label_subset,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
print(plot)

celltype_list = sort(unique(as.character(Idents(data_harmony_run_label_subset))))
celltype_list

celltype_num = sort(table(as.character(Idents(data_harmony_run_label_subset))))

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','Plasma Cell')

celltype_list=  c('Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','Plasma Cell')
celltype_list = c('CD14+ Mono')
celltype = celltype_list[1]


DE_type = 'DESeq2'
DE_type = 'EdgeR'

data_harmony_run_label_subset$Sample_type[DE_input$'Sample Type' == "PBMC"] = 'PB'
data_harmony_run_label_subset$Sample_type[DE_input$'Sample Type' == "Bone Marrow"] = 'BM'


for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0('baseline BM ',celltype)
  ident2 = paste0('baseline PB ',celltype)
  
  DE_input = data_harmony_run_label_subset
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Sample_type = ''
  DE_input$Sample_type[DE_input$'Sample Type' == "PBMC"] = 'PB'
  DE_input$Sample_type[DE_input$'Sample Type' == "Bone Marrow"] = 'BM'
  

  
  DE_input$DE_ident = paste0(DE_input$Treatment,' ', DE_input$Sample_type, ' ', Idents(DE_input))


  DE_input = DE_input[,Idents(DE_input) %in% c(celltype)]
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
  ncol(DE_input)
  unique(DE_input$DE_ident)
  

  IdentPerPatient = table(DE_input$`Patient Number`,DE_input$DE_ident )
  print(IdentPerPatient)
  #patient_names = rownames(IdentPerPatient)
  
  
  #patient_keep = patient_names[rowSums(IdentPerPatient >= 30) == 2]
  #patient_keep = patient_names[rowSums(IdentPerPatient >= 30) >= 1]
  
  IdentPerSample = table(DE_input$sample )
  
  sample_list = names(IdentPerSample)
  sample_keep = sample_list[IdentPerSample >= 20]
  # 

    # if (grepl('NBM', ident1, fixed = TRUE) | grepl('NBM', ident2, fixed = TRUE) |
  #     grepl(' PR ', ident2, fixed = TRUE)  ){
  #   patient_keep = patient_names[rowSums(IdentPerPatient >= 30) == 1]
  # }else{
  #   patient_keep = patient_names[rowSums(IdentPerPatient >= 30) == 2]
  # }
  if (length(sample_keep) < 2){
    print('Not enough samples')
    next
  }
  DE_input = DE_input[,DE_input$sample %in% sample_keep]
  if (ncol(DE_input) < 100){
    print('Not enough cells')
    next
  }
  
  print(sample_keep)
  
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
    print('Idents not present')
    next
  }

  if (length(unique(kit)) == 1){
    formula = ~DE_ident + DR + Patient
  }else{
    formula = ~DE_ident + DR + kit + Patient
  }
  
  
  design <- model.matrix(formula)
  colnames(design) <- gsub("DE_ident", "", colnames(design))
  colnames(design)
  is.fullrank(design)
  design_fr = fullRank(design)
  is.fullrank(design_fr)
  colnames_old = colnames(design)
  colnames_fr = colnames(design_fr)
  
  #print(colnames_old[!(colnames_old %in% colnames_fr)])
  # Make contrast matrix that says you want to compare ident2 to ident1
  int2 = match(ident2,colnames(design_fr) )
  con <- integer(ncol(design_fr))
  con[int2] <- 1
  
  # Filter out poorly expressed genes
  keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                       min.count = 1,min.total.count=10, 
                       large.n = 10,min.prop = 0.1)
  print(sum(keep))
  #next
  
  subfolder = paste0(ident1,' Vs ',ident2)
  
  if (DE_type == 'DESeq2'){
    result_DESeq2 = runDESeq2(DE_input,design_fr,contrast = con, keep = keep, 
                              folder_output = filepath_cluster, subfolder = subfolder)
    
    
    if (all(is.na(result_DESeq2))){
      next
      }
    dds = result_DESeq2[[1]]
    res = result_DESeq2[[2]]
    res$logFC = res$log2FoldChange
  }else if ( DE_type == 'EdgeR'){
  
    result_edgeR = runEdgeR(DE_input,design_fr, contrast = con, keep = keep,
                            folder_output = filepath_cluster,
                            subfolder = subfolder)
    
    dge_edgeR = result_edgeR[[1]]
    fit_edgeR = result_edgeR[[2]]
    qlf = result_edgeR[[3]]
    res = result_edgeR[[4]]
    res$pvalue = res$PValue
    res$log2FoldChange = res$logFC
    res$gene = rownames(res)
  }
  
  names(res)[names(res) == "logFC"] <- "log2FoldChange"
  names(res)[names(res) == "FDR"] <- "padj"
  names(res)[names(res) == "PValue"] <- "pvalue"
  
  
  
  res = res[order(res$log2FoldChange),]
  res = res[res$pvalue < 0.05,]
  
  
  folder = paste0(filepath_cluster,'/DE/',DE_type,'/',subfolder,'/')
  dir.create(folder,recursive = T)
  
  filename = paste0(folder,'/DE_',DE_type,' ',subfolder,'_clean','.csv')
  
  write.csv(res, file = filename,row.names=TRUE)
  
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl("^RP[SL]", res$gene),]
  

  max_gene = 50
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
  

  
}
########################
## Per Patient DE
########################

patientByGroup =  table(data_harmony_run_label_subset$`Patient Number`,
                        paste0(data_harmony_run_label_subset$Treatment,' ', data_harmony_run_label_subset$'Sample Type'))

patientByGroup = as.data.frame(patientByGroup > 0)
patient_list = rownames(patientByGroup)[patientByGroup$`baseline Bone Marrow` & patientByGroup$`baseline PBMC`]

#patient_list = sort(unique(data_harmony_run_label$`Patient Number`))
#patient_list

patient = patient_list[1]

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'dMono','Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype_list=  c('Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

#data_harmony_run_label = renameCells(data_harmony_run_label,idents = c('cDC1','cDC2'),newident = 'DC')
for (patient in patient_list){
  print(patient)
  for (celltype in celltype_list){
    print(celltype)
    #celltype = 'eTreg'
    ident1 = paste0('baseline BM ',celltype)
    ident2 = paste0('baseline PB ',celltype)
    
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
    
    data = as.data.frame(DE_input@assays[["RNA"]]@counts)
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
      print('Not all idents present')
      next
    }
    library(robustbase)
    library(DESeq2)
    formula = ~  DE_ident + DR
    
    #formula = ~  DE_ident + DR + kit
    design <- model.matrix(formula)
    colnames(design) <- gsub("DE_ident", "", colnames(design))
    colnames(design)
    is.fullrank(design)
    design_fr = fullRank(design)
    is.fullrank(design_fr)
    colnames_old = colnames(design)
    colnames_fr = colnames(design_fr)
    
    print(colnames_old[!(colnames_old %in% colnames_fr)])
    
    #int1 = match(ident1,colnames(design_fr) )
    int2 = match(ident2,colnames(design_fr) )
    con <- integer(ncol(design_fr))
    #con[int1] <- 1 
    con[int2] <- 1
    
    
    
    keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                         min.count = 1,min.total.count=10, 
                         large.n = 10,min.prop = 0.1)
    print(sum(keep))

    subfolder = paste0('Patient',patient,'/',ident1,' Vs ',ident2)
    
    
    colnames_old[!(colnames_old %in% colnames_fr)]
    

    result_edgeR = runEdgeR(DE_input,design_fr, contrast = con, keep = keep,
                            folder_output = filepath_cluster,
                            subfolder = subfolder)
    
    dge_edgeR = result_edgeR[[1]]
    fit_edgeR = result_edgeR[[2]]
    qlf = result_edgeR[[3]]
    res = result_edgeR[[4]]
    res$pvalue = res$PValue
    res$log2FoldChange = res$logFC
    res$gene = rownames(res)
    
    res = res[order(res$log2FoldChange),]
    res = res[res$pvalue < 0.05,]
    
    res = res[!grepl("MT-", res$gene),]
    res = res[!grepl(" ?RP\\w+ ?", res$gene),]
    
    folder = paste0(filepath_cluster,'/DE/EdgeR/',subfolder,'/')
    dir.create(folder,recursive = T)
    filename = paste0(folder,'Patient',patient,'_clean','.csv')
    
    write.csv(res, file = filename,row.names=TRUE)
    
    #png(file=pathName, height=12, width=20)
    #png(file=pathName)
    max_gene = 50
    res$gene = as.character(res$gene)
    res_pos = res[res$log2FoldChange > 0,]
    res_neg = res[res$log2FoldChange < 0,]
    
    res_pos = res_pos[order(-res_pos$logFC),]
    
    
    if (nrow(res_pos) > max_gene){
      
      res_pos = res_pos[1:max_gene,]
      res_pos = res_pos[rowSums(is.na(res_pos)) != ncol(res_pos),]
      
    }
    if (nrow(res_neg) > max_gene){
      
      res_neg = res_neg[1:max_gene,]
      res_neg = res_neg[rowSums(is.na(res_neg)) != ncol(res_neg),]
      
    }


    
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
  }
}
