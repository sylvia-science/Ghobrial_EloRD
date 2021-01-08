source('~/Desktop/scRNA/Code/Analysis/DE_Methods.R') # Chance to where your DE_Methods file is
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

for (celltype in celltype_list){
  print(celltype)
  ident1 = paste0('EOT GR ',celltype)
  ident2 = paste0('EOT PR ',celltype)
  
  DE_input = data_merge_run_label
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response == 'MR' ] = 'PR'
  DE_input$Best_Overall_Response[DE_input$Best_Overall_Response %in% c('VGPR','CR','sCR') ] = 'GR'
  
  DE_input$DE_ident = paste0(DE_input$Treatment, ' ', 
                             DE_input$Best_Overall_Response, ' ', Idents(DE_input))
  DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)] #Subset to only the idents you're looking at
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
  
  if (length(unique(DE_ident)) == 1){ # If not all idents exist, don't continue
    next
  }
  

  formula = ~DE_ident + DR + kit + Patient
  
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
  
  # This is where you tell it what ident you're interested in. 
  # The design will only have one of your idents available to compare against the other
  # This ident could by ident1 or ident 2
  #int1 = match(ident1,colnames(design_fr) ) 
  int2 = match(ident2,colnames(design_fr) )
  con <- integer(ncol(design_fr))
  #con[int1] <- 1 
  con[int2] <- 1
  
  
  
  # Filtering for more expressed genes
  keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                       min.count = 1,min.total.count=10, 
                       large.n = 10,min.prop = 0.1)
  print(sum(keep))
  #next
  
  subfolder = paste0(ident1,' Vs ',ident2)
  
  result_DESeq2 = runDESeq2(DE_input,design_fr,contrast = con, keep = keep, 
                              folder_output = filepath_cluster, subfolder = subfolder)
    
    
  if (all(is.na(result_DESeq2))){
    next
  }
  dds = result_DESeq2[[1]]
  res = result_DESeq2[[2]]
  res$logFC = res$log2FoldChange

  res = res[order(res$log2FoldChange),]
  res = res[res$pvalue < 0.05,]
  
  # Remove MT and RP genes
  res = res[!grepl("MT-", res$gene),]
  res = res[!grepl(" ?RP\\w+ ?", res$gene),]
  
  # Make a heatmap of the DE genes
  folder = paste0(filepath_cluster,'/DE/',DE_type,'/',subfolder,'/')
  dir.create(folder,recursive = T)
  pathName <- paste0(folder,paste0('HeatMap','.png'))

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