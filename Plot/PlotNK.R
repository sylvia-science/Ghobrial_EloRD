X = as.data.frame(t(data_matrix))

umap_new = as.data.frame(output)
umap_new$cluster = umap_new$label

clusters_to_keep = c("Tgd"  ,"aCX3CR1+ GZMB+ CD56dim" ,   "aCXCR4+ CD56dim"    
                                 ,"aGZMK+ CCL3+ CD56dim",  "cCD56dim" ,"CD56bright"   ,"NFkB-high",'aCCL3+ GZMK+')  
markers <- c("CD3D","CD3G","NKG7","TRGC1","TRGC2","NCAM1","SELL","CD2","GZMK","GZMB","GZMH","FCGR3A","CX3CR1","PRF1","CCL3","JUN","FOSB","FOS","DUSP1","JUNB","IER2","CXCR4","REL","RELB","NFKB1","NFKB2","BIRC2","AREG","KDM6B","TNFAIP3","NFKBIA","TGFB1","FAM177A1","CD160","XCL1","XCL2","CD52","CD122","CD27","CD69","LTB","CD44","CD127","GATA3","RORC","NCR1", "CD4","LY6A","KLRB1C","RORA","CD83","MAFF", "IL2RA","CCR7")
marker_counts <- X[,colnames(X) %in% markers]
marker_counts$idents <- umap_new$ident[match(rownames(marker_counts), rownames(umap_new))]
marker_means <- aggregate(marker_counts[, 1:(ncol(marker_counts)-1)], list(marker_counts$idents), mean)
rownames(marker_means) <- paste0("C",marker_means$Group.1)

marker_means$Group.1  <- NULL
marker_means <- data.frame(t(marker_means))
#marker_means <- marker_means[,colnames(marker_means) %in% paste0("C",clusters_to_keep)]


marker_means$Marker <- rownames(marker_means)
marker_means$cyto <- rowMeans(marker_means[,c("C2","C8","C16")])
marker_means$nfkb <- rowMeans(marker_means[,c("C5","C14")])
marker_means$C2 <- NULL
marker_means$C8 <- NULL
marker_means$C16 <- NULL
marker_means$C17 <- NULL
marker_means$C19 <- NULL
marker_means$C5 <- NULL
marker_means$C14 <- NULL
#marker_means$Diff <- marker_means$C1 - marker_means$C0
marker_means_long <- gather(marker_means, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = c("C3","C4","cyto","C6","C9","C11","nfkb"))
marker_means_long = marker_means_long[!is.na(marker_means_long$Cluster),]
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
marker_means_long$Cluster <- as.character(marker_means_long$Cluster)
marker_means_long$Cluster[marker_means_long$Cluster == "C4"] <- "CD56br"
marker_means_long$Cluster[marker_means_long$Cluster == "C3"] <- "Tgd"
marker_means_long$Cluster[marker_means_long$Cluster == "cyto"] <- "cCD56dim"
marker_means_long$Cluster[marker_means_long$Cluster == "C6"] <- "aCX3CR1+GZMB+CD56dim"
marker_means_long$Cluster[marker_means_long$Cluster == "C11"] <- "aGZMK+CCL3+CD56dim"
marker_means_long$Cluster[marker_means_long$Cluster == "C9"] <- "aCXCR4+CD56dim"
marker_means_long$Cluster[marker_means_long$Cluster == "nfkb"] <- "NFkB-high"
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels = c("Tgd","CD56br", "cCD56dim","aGZMK+CCL3+CD56dim", "aCXCR4+CD56dim","aCX3CR1+GZMB+CD56dim","NFkB-high"))

pdf(paste0(filepath_cluster,'Plots/NK_Annotation_Heatmap.pdf'), width=12,height = 4)
plot = ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + scale_fill_gradient(low="white",high="orangered3") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Exp") +theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) + ylab("")
print(plot)
dev.off()

