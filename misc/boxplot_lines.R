
data = data_run_subset_label
sample_list = unique(data$sample)
celltype_list = unique(Idents(data))
stats_summary = data.frame(matrix(ncol = length(celltype_list), nrow = length(sample_list)))
colnames(stats_summary) = celltype_list
rownames(stats_summary) = sample_list

stats_summary = (table(data_run_subset_label$sample,Idents(data_run_subset_label)))
stats_summary = stats_summary/rowSums(stats_summary)

cell = 'Th17'
stats_summary$Sample = rownames(stats_summary)
stats_summary =merge(stats_summary, metaData, by = 'Sample')
stats_summary$Response = metaData$Best_Overall_Response[metaData$Sample %in% stats_summary$Sample ]
stats_summary$category = as.character(stats_summary$Treatment)

category_levels = c('baseline','C9D1','EOT' )
stats_summary_line = stats_summary[stats_summary$category %in% category_levels,]    
stats_summary_line$category = factor(stats_summary_line$category, levels  = category_levels)
stats_summary_line$Patient = stats_summary_line$`Patient Number`


plot = ggplot(stats_summary_line, aes(x = category, y = !!ensym(cell))) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, max_val))+
  xlab("") + ylab(paste0(cell, ' Proportion'))+
  geom_point(color="black", size=2) +
  geom_line(aes(group=Patient, color=stats_summary_line$Response,alpha = Response ), size = 1) +
  #scale_colour_manual(values=c(category1="blue",category2="red"))+
  scale_alpha_manual(values = c( 0.9, 0.9))+
  guides(alpha = FALSE)+
  labs(color = "Response")+
  theme_classic()

plot = plot + theme(
  axis.title.x = element_text(color="black", size=24 ),
  axis.title.y = element_text(color="black", size=24),
  axis.text= element_text(color="black", size=24),
  legend.text=element_text(size=18),
  legend.title=element_text(size=18)
  
)

print(plot)


  
  