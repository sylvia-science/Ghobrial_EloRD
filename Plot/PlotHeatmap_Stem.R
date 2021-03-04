#Plot heatmap of R markers 
library(tidyr)
celltype = c('HSC','MDPC', 'GMPC')

celltype = c('HSC','CMPC','MDPC', 'GMPC','Pro B-cell', 'Pre B-cell')

unique((Idents(data_harmony_run_label)))
data_input  =data_harmony_run_label[ ,Idents(data_harmony_run_label) %in% celltype]

X = as.data.frame(t(data_input@assays[["RNA"]]@counts))


markers = c('CDK6','CD34','SPINK2', 'TSC22D1', #HSC
            'CD34', 'CDK6', 'FLT3', 'KIT', 'MIF', 'NPM1', #CMPC
            'STMN1', 'TUBB', 'BIRC5', #MDPC
            'MPO', 'PRTN3', 'ELANE', 'AZU1', 'CTSG', 'LYZ', #GMPC
            'CD79B', 'DNTT', 'VPREB1', 'CD9', #Pro B
            'DUT', 'PCNA', 'IGLL1','CD24') #Pre B


markers = c('KIAA0125','PRSS57','CD34', 'CDK6', #HSC
             'FLT3', 'KIT', 'MIF', 'NPM1', #CMPC
             'TOP2A', 'MKI67', 'PCNA','TUBB','LYZ', 'SRGN', 'FCN1', 'S100A4','KLF4', 'NR4A1', 'CSF1R', 'FCGR3A', 'FCGR2A', #MDPC
             'ELANE', 'PRTN3', 'MPO', 'AZU1', #GMPC
            'CD19', 'MS4A1', 'IGHM', 'MME', 'DNTT', 'CD79A', 'CD37', 'CD24', 'IGLL1', 'VPREB1' #Pre/Pro B
            ) 

markers = c('CD34', 'CDK6', 'FLT3', 'KIT', 'MIF', 'NPM1' ,'LYZ', 'FCN1', 'S100A4', 'KLF4', 'CSF1R', 'FCGR2A', 'SRGN', 
            'ELANE', 'PRTN3', 'MPO', 'AZU1', 'PRSS57', 'MME', 'DNTT', 
            'CD19', 'MS4A1', 'IGHM', 'CD79A', 'CD24', 'IGLL1')
umap_new = as.data.frame(Idents(data_input))
colnames(umap_new) = 'idents'
marker_counts <- X[,colnames(X) %in% markers]

marker_counts$idents <- umap_new$idents[match(rownames(marker_counts), rownames(umap_new))]
#marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
marker_means <- aggregate(scale(marker_counts[, 1:(ncol(marker_counts)-1)]), list(marker_counts$idents), mean)

#rownames(marker_means) <- paste0("C",marker_means$Group.1)

colname = marker_means$Group.1 
marker_means$Group.1  <- NULL
marker_means <- data.frame(t(marker_means))
colnames(marker_means) =  colname
marker_means$Marker = rownames(marker_means)


#marker_means$Diff <- marker_means$C1 - marker_means$C0
marker_means_long <- gather(marker_means, "Cluster", "Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = celltype)
marker_means_long = marker_means_long[!is.na(marker_means_long$Cluster),]
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)


pdf(paste0(filepath_cluster,'Plots/heatmap_StemCells.pdf'), width=12,height = 4)
plot = ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Exp)) +
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Exp),0,max(marker_means_long$Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) + ylab("")
print(plot)
dev.off()


