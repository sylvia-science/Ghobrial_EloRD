prop = read.table('/home/sujwary/Desktop/Final_proportion_table_normalized.txt',check.names = F)

old_colnames = sort(colnames(prop))
old_colnames = gsub('_',' ' ,old_colnames)
old_colnames = gsub("[[:punct:]]",'' ,old_colnames)
old_colnames


remove_list = c('Remove','14','32','Erythrocyte','DC/T-Cell DBL',
                'Mono/CD8+ T Cell DBL','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL','42','41','',0:40
                ,'CD14+ Mono/T-cell DBL','CD14+ Mono/CD8+ T-cell DBL','Pro Erythrocyte',
                'dDC','DC/T-cell DBL','dMIP1a+ Mono','sMono','MK','dMono','sCD14+ Mono','Remove 29',
                'Remove 30','dT-cell','CD16+ Mono/T-cell DBL')

data_harmony_run_label_remove = data_harmony_run_label[,!(Idents(data_harmony_run_label) %in% remove_list)]
Idents(data_harmony_run_label_remove) = as.character(Idents(data_harmony_run_label_remove))


new_colnames = sort(as.character(sort(unique(Idents(data_harmony_run_label_remove)))))
new_colnames = gsub("[[:punct:]]",'' ,new_colnames)
new_colnames = gsub("Tcell",'T cell' ,new_colnames)


old_colnames[!(old_colnames %in% new_colnames)]

new_colnames[!(new_colnames %in% old_colnames)]


length(old_colnames)
length(new_colnames)

table = (table( data_harmony_run_label_remove$sample,Idents(data_harmony_run_label_remove)))

write.csv(table, file = paste0(filepath_cluster,'/data/','sampleByIdent_Table.csv'))

