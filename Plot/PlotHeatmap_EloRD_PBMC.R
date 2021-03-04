
library(pheatmap)
library(tidyr)

#Plot heatmap of NK Cells

celltype = c("Tgd","CD56bright","aCCL3+ CD56dim","cCD56dim","NFkB-high")
celltype_factor = c("Tgd","CD56bright NK","cCD56dim NK","aCCL3+ GZMK+ NK","NFkB-high NK")
celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]

data_input  =data_run_subset_label[ ,Idents(data_run_subset_label) %in% celltype]
unique(Idents(data_run_subset_label))
unique(Idents(data_input))
Ident_list = as.character(Idents(data_input))

Ident_list[Ident_list == "CD56bright"] <- "CD56bright NK"
Ident_list[Ident_list == "Tgd"] <- "Tgd"
Ident_list[Ident_list == "cCD56dim"] <- "cCD56dim NK"
Ident_list[Ident_list== "NFkB-high"] <- "NFkB-high NK"
Ident_list[Ident_list == "aCCL3+ CD56dim"] <- "aCCL3+ GZMK+ NK"
unique(Ident_list)

mat = as.data.frame(t(data_input@assays[["RNA"]]@counts))

markers <- c("CD3D","CD3G","TRGC1","TRGC2","NCAM1","SELL","CD2","GZMK","GZMB","GZMH","FCGR3A","CX3CR1","PRF1","CCL3","JUN","FOSB","FOS","DUSP1","JUNB","IER2","CXCR4","REL","RELB","NFKB1","NFKB2","BIRC2","AREG","KDM6B","TNFAIP3","NFKBIA","TGFB1")
marker_counts <- mat[,colnames(mat) %in% markers]
rownames(marker_counts) <- gsub("\\.","-",rownames(marker_counts))
marker_counts$ident <- Ident_list

marker_means <- aggregate(scale(marker_counts[, 1:(ncol(marker_counts)-1)]), list(marker_counts$ident), mean)
rownames(marker_means) = marker_means$Group.1
marker_means$Group.1 = NULL
marker_means <- as.data.frame(t(marker_means))
marker_means$Marker <- rownames(marker_means)
marker_means_long <- gather(marker_means, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- as.character(marker_means_long$Cluster)
unique(marker_means_long$Cluster )

marker_means_long$Cluster <- factor(marker_means_long$Cluster,levels = celltype_factor)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)

pdf("/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Plots/Paper/Blood_NK heatmap.pdf",width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Exp") +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()

########################################################################
# T cells
celltype = c("Naive CD4+ T-cell","Stim Naive CD4+ T-cell","IFN+ CD4+ T-cell","TSCM","CD4+ TCM",
             "cTreg","eTreg","Th2","TRM","Th17","aTh17",
             "Naive CD8+ T-cell","CD8+ TCM","GZMK+ CD8+ T-cell","GZMK+ CCL3+ CCL4+ CD8+ T-cell","GZMH+ GZMB+ CD8+ T-cell","TEMRA")
celltype_factor = c("Naive_CD4_T_cell","Stim_Naive_CD4_T_cell","IFN_CD4_T_cell","TSCM","CD4_TCM",
                    "cTreg","eTreg","Th2","TRM","Th17","aTh17",
                    "Naive_CD8_T_cell","CD8_TCM","GZMK_CD8_T_cell","GZMK_CCL3_CCL4_CD8_T_cell","GZMH_GZMB_CD8_T_cell","TEMRA")

data_input  =data_run_subset_label[ ,Idents(data_run_subset_label) %in% celltype]
unique(Idents(data_run_subset_label))
unique(Idents(data_input))
Ident_list = as.character(Idents(data_input))
celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]
#unique(Idents(data_input))[!( unique(Idents(data_input)) %in% celltype)]
mat = as.data.frame(t(data_input@assays[["RNA"]]@counts))


markers <- c("CD4","SELL","CCR7","JUN","FOSB","DUSP1","ISG15", "IFI6","IFI44L","IL7R","BCL2","CDK6","FOXP3","CTLA4","HLA-DRB1","CD74","GATA3","ITGB1","CRIP1","LGALS1","S100A10","S100A11","RORA","KLRB1","CCL5","GZMK","CD8A","CD8B","CCL3","CCL4","GZMH","GZMB","PRF1","GNLY","FCGR3A")
marker_counts <- mat[,colnames(mat) %in% markers]
rownames(marker_counts) <- gsub("\\.","-",rownames(marker_counts))
marker_counts$ident <- Ident_list
marker_counts$ident <- gsub("\\+","",marker_counts$ident )
marker_counts$ident <- gsub(" ","_",marker_counts$ident )
marker_counts$ident <- gsub("-","_",marker_counts$ident )
marker_means <- aggregate(scale(marker_counts[, 1:(ncol(marker_counts)-1)]), list(marker_counts$ident), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means <- data.frame(t(marker_means))
marker_means$Marker <- rownames(marker_means)
marker_means_long <- gather(marker_means, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels=celltype_factor)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)

pdf("/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Plots/Paper/Blood_T heatmap.pdf",width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()
###############################################################
## Monocytes

celltype = c("SELL+ CD14+ Mono","TGFb1+ CD14+ Mono","CD14+ CD16+ Mono","CD16+ Mono","IFN+ Mono")
celltype_factor = c("SELL+ CD14+ Monocytes","TGFb1+ CD14+ Monocytes","CD14+ CD16+ Monocytes","CD16+ Monocytes","IFN+ Monocytes")

data_input  =data_run_subset_label[ ,Idents(data_run_subset_label) %in% celltype]
unique(Idents(data_run_subset_label))
unique(Idents(data_input))
Ident_list = as.character(Idents(data_input))
Ident_list[Ident_list == "SELL+ CD14+ Mono"] <- "SELL+ CD14+ Monocytes"
Ident_list[Ident_list == "TGFb1+ CD14+ Mono"] <- "TGFb1+ CD14+ Monocytes"
Ident_list[Ident_list == "CD14+ CD16+ Mono"] <- "CD14+ CD16+ Monocytes"
Ident_list[Ident_list == "CD16+ Mono"] <- "CD16+ Monocytes"
Ident_list[Ident_list == "IFN+ Mono"] <- "IFN+ Monocytes"

celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]
mat = as.data.frame(t(data_input@assays[["RNA"]]@counts))


markers <- c("CD14","SELL","S100A8", "S100A9", "S100A12","HLA-DRA","HLA-DRB1","TIMP1","TGFB1","CD300E","LGALS3","FCGR3A","CSF1R","MS4A7","CD74","MX1","ISG15","IFI6","IFI44L")
marker_counts <- mat[,colnames(mat) %in% markers]
marker_counts$ident <- Ident_list
marker_means <- aggregate(scale(marker_counts[, 1:(ncol(marker_counts)-1)]), list(marker_counts$ident), mean)
rownames(marker_means) = marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means <- as.data.frame(t(marker_means))
marker_means$Marker <- rownames(marker_means)
marker_means_long <- gather(marker_means, "Cluster", "Mean_Exp",-Marker)

marker_means_long$Cluster <- factor(marker_means_long$Cluster,levels = celltype_factor)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)

pdf("/disk2/Projects/EloRD/Output/Harmony/AllSamples_PBMC/Batch_Sample_Kit/Plots/Paper/Blood_Mono heatmap.pdf",width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp))))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Exp") +theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()


