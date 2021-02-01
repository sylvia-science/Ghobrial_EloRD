
data_input = data_run_subset
cell_names_all = sub("_.*", "", colnames(data_input))

data_input$sample_cell = paste0(data_input$Sample,' ',cell_names_all)

path = paste0(filepath_cluster,'/Data/')
dir.create(path)
labels = as.character(Idents(data_input))


labels = as.data.frame(colnames(data_input))
colnames(labels) = 'Cell'
labels$sample = (as.character(data_input$sample))

labels$cell_type = (as.character(Idents(data_input)))
labels$sample_cell = data_input$sample_cell

write.csv(labels, file = paste0(path,'Labels.csv'), )

