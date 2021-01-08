source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R')
library(robustbase)
library(DESeq2)
library(EnrichmentBrowser)



data_merge_run_label = ScaleData(data_merge_run_label, features = rownames(data_merge_run_label))
data_merge_run_label = renameCells(data_merge_run_label,idents = c('cDC1','cDC2'),newident = 'DC')
data_merge_run_label = renameCells(data_merge_run_label,idents = c('CD14+CD16+ Mono'),newident = 'CD16+ Mono')

celltype_list = sort(unique(as.character(Idents(data_merge_run_label))))
celltype_list

celltype_num = sort(table(as.character(Idents(data_merge_run_label))))

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
celltype_list = c('aTh17')
celltype = celltype_list[1]


DE_type = 'DESeq2'
DE_type = 'EdgeR'

for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0('EOT GR ',celltype)
  ident2 = paste0('EOT PR ',celltype)
  
  #ident1 = paste0('CD16+ Mono')
  #ident2 = paste0('CD14+ Mono')
  
  DE_input = data_merge_run_label
  #DE_input = renameCells(DE_input,idents = c('TIMP1+ CD14+ Mono','SELL+ CD14+ Mono'),
  #                       newident = 'CD14+ Mono')
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
                            DE_input$Best_Overall_Response, ' ', Idents(DE_input))
  #DE_input$DE_ident = paste0(DE_input$Treatment,' ',DE_input$'Dexa or not', ' ', Idents(DE_input))
  #DE_input$DE_ident = paste0(DE_input$Treatment,' ', Idents(DE_input))
  
  #DE_input$DE_ident = paste0(Idents(DE_input))
  #DE_input = DE_input[,DE_input$Treatment =='baseline']
  DE_input = DE_input[,Idents(DE_input) %in% c(celltype)]
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
  ncol(DE_input)
  unique(DE_input$DE_ident)
  

  IdentPerPatient = table(DE_input$`Patient Number`,DE_input$DE_ident )
  print(IdentPerPatient)
  patient_names = rownames(IdentPerPatient)
  
  if (grepl('NBM', ident1, fixed = TRUE) | grepl('NBM', ident2, fixed = TRUE) |
      grepl(' PR ', ident2, fixed = TRUE)  ){
    patient_keep = patient_names[rowSums(IdentPerPatient >= 30) == 1]
  }else{
    patient_keep = patient_names[rowSums(IdentPerPatient >= 30) == 2]
  }
  if (length(patient_keep) < 2){
    print('Not enough patients')
    next
  }
  DE_input = DE_input[,DE_input$`Patient Number` %in% patient_keep]
  if (ncol(DE_input) < 100){
    print('Not enough cells')
    next
  }
  
  print(patient_keep)
  
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

  if (length(unique(kit)) == 1){
    formula = ~DE_ident + DR + Patient
  }else{
    formula = ~DE_ident + DR + kit + Patient
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
  res = res[order(res$log2FoldChange),]
  res = res[res$pvalue < 0.05,]
  
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl("^RP[SL]", res$gene),]
  
  
  folder = paste0(filepath_cluster,'/DE/',DE_type,'/',subfolder,'/')
  dir.create(folder,recursive = T)
  pathName <- paste0(folder,paste0('HeatMap','.png'))
  #png(file=pathName, height=12, width=20)
  #png(file=pathName)
  len = nrow(res)
  if (nrow(res) > 200){
    res = res[c(1:100,(len - 100):len),]
    res = res[rowSums(is.na(res)) != ncol(res),]
    res = unique(res)
  }
  
  
  png(file=pathName,width=2000, height=40*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input, features = res$gene, group.by = 'DE_ident')
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
## Per Patient DE
########################

patient_list = sort(unique(data_merge_run_label$`Patient Number`))
patient_list

patient = patient_list[14]

celltype_list=  c('CD14+ Mono','TRM','Naive CD4+ T-cell','Naive CD8+ T-cell','GZMK+ CD8+ T-cell',
                  'cTreg', 'Mature NK','CD56Br NK','Th17','Inhibitory NK','TEMRA','DC',
                  'dMono','Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

celltype_list=  c('Intermediate CD4+ T-cell','CD4+ TCM','B Cell','GZMK+ CCL3+ CCL4+ CD8+ T-cell',
                  'GZMH+ GZMB+ CD8+ T-Cell','Th2','SELL+ CD14+ Mono','CCL5+ CD4+ T-cell','aTh17',
                  'CD16+ Mono','sMono','Cytotoxic NK','ILC1')

data_merge_run_label = renameCells(data_merge_run_label,idents = c('cDC1','cDC2'),newident = 'DC')

for (celltype in celltype_list){
  print(celltype)
  #celltype = 'eTreg'
  ident1 = paste0('baseline ',celltype)
  ident2 = paste0('C9D1 ',celltype)
  
  #ident1 = paste0('CD16+ Mono')
  #ident2 = paste0('CD14+ Mono')
  
  DE_input = data_merge_run_label
  DE_input = DE_input[,DE_input$`Patient Number` == patient]
  
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  
  DE_input$DE_ident = paste0(DE_input$Treatment,' ', Idents(DE_input))
  
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
  
  colnames_old[!(colnames_old %in% colnames_fr)]
  
  keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                       min.count = 1,min.total.count=10, 
                       large.n = 10,min.prop = 0.1)
  subfolder = paste0('Patient',patient, '/', ident1,' Vs ',ident2)
  
  result_DESeq2 = runDESeq2(DE_input,design_fr,contrast, keep, 
                            folder_output = filepath_cluster, subfolder = subfolder)
  
  dds = result_DESeq2[[1]]
  res = result_DESeq2[[2]]
  #plotDispEsts(dds)
  
  res = res[order(res$log2FoldChange),]
  res = res[res$pvalue < 0.05,]
  
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl(" ?RP\\w+ ?", res$gene),]
  
  folder = paste0(filepath_cluster,'/DE/DESeq2/',subfolder,'/')
  dir.create(folder,recursive = T)
  pathName <- paste0(folder,paste0('HeatMap','.png'))
  #png(file=pathName, height=12, width=20)
  #png(file=pathName)
  len = nrow(res)
  if (len > 200){
    res = res[c(1:100,(len - 100):len),]
    res = res[rowSums(is.na(res)) != ncol(res),]
    res = unique(res)
  }
  
  print(nrow(res))
  png(file=pathName,width=2000, height=40*nrow(res),res = 100)
  
  plot = DoHeatmap(DE_input, features = res$gene, group.by = 'DE_ident')
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