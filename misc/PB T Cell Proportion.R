# Get mean num of PB T cells and mean proportion

# Remove junk cells and take only baseline PBMCs
data_harmony_run_label_remove = data_harmony_run_label[,data_harmony_run_label$`Sample Type` == 'PBMC']
data_harmony_run_label_remove = data_harmony_run_label_remove[,!(Idents(data_harmony_run_label_remove) %in% c('Remove', '13','30'))]
data_harmony_run_label_remove = data_harmony_run_label_remove[,data_harmony_run_label_remove$Treatment == 'baseline']
Idents(data_harmony_run_label_remove) = as.character(Idents(data_harmony_run_label_remove))

# Add CD8 T cells to CD4 T cells
idents = as.character(Idents(data_harmony_run_label_remove))
idents[idents == 'CD8+ T-cell'] = 'T-cell'

Idents(data_harmony_run_label_remove)  = idents
#plot = DimPlot(data_harmony_run_label_remove,pt.size = 0.7, reduction = "umap",label = TRUE,label.size = 8)
#print(plot)


num_Cells = table(Idents(data_harmony_run_label_remove), data_harmony_run_label_remove$sample)

write.csv(as.data.frame.matrix(num_Cells), '/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Cluster/PCA40/res3/Data/num_Cells.csv')

num_Cells

mean_num_T_Cell = mean(num_Cells['T-cell',])

proportions =  num_Cells['T-cell',]/colSums(num_Cells)

colSums(num_Cells)

mean_proportion_T_Cell = mean(proportions)

