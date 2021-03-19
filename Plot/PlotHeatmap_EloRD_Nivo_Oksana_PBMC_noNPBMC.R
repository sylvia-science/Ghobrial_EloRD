
library(pheatmap)
library(tidyr)
library("readr")

#Plot heatmap of NK Cells

NFkB = as.character(unique(Idents(data_run_subset_label))[11])
celltype =c ("CD56br NK","aCD56br NK",NFkB,"CD56dim NK","aCD56dim NK","CCL3+ NK","IFN+ NK")
celltype_factor = c("CD56br NK","aCD56br NK",NFkB,"CD56dim NK","aCD56dim NK","CCL3+ NK","IFN+ NK")
celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]
unique(Idents(data_run_subset_label))[!(unique(Idents(data_run_subset_label)) %in% celltype)]

unique(Idents(data_run_subset_label))


data_input  =data_run_subset_label[ ,Idents(data_run_subset_label) %in% celltype]
unique(Idents(data_input))
Ident_list = as.character(Idents(data_input))

Ident_list[Ident_list == "CD56bright"] <- "CD56bright NK"
Ident_list[Ident_list == "Tgd"] <- "Tgd"
Ident_list[Ident_list == "cCD56dim"] <- "cCD56dim NK"
Ident_list[Ident_list== "NFkB-high"] <- "NFkB-high NK"
Ident_list[Ident_list == "aCCL3+ CD56dim"] <- "aCCL3+ GZMK+ NK"
Ident_list[Ident_list == "CD81+ cCD56dim"] <- "CD81+ cCD56dim NK"

unique(Ident_list)

mat = as.data.frame(t(data_input@assays[["RNA"]]@data))

markers <- c("SELL","CD2","NCAM1","XCL2","KLRC1","JUN","FOSB","DUSP1",
             "CXCR4","IFRD1","AREG","NFKB1","NFKB2","REL","RELB","KDM6B","TGFB1","TNFAIP3","NFKBIA",
             "GZMB","GZMH","FCGR3A","PRF1","CX3CR1","CCL3","CCL3L3","CCL4","CCL4L2",
             "ISG15","IFI6","IFI44L")
marker_counts <- mat[,colnames(mat) %in% markers]
rownames(marker_counts) <- gsub("\\.","-",rownames(marker_counts))
marker_counts$ident <- Ident_list

colnames(marker_counts)[1:3] 

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

pdf(paste0('/disk2/Projects/EloRD/Paper/EloRD_Nivo_Oksana/',folder_name,'_NK heatmap.pdf'),width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),
                       values=scales::rescale(c(min(marker_means_long$Mean_Exp), 0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()


########################################################################
# T cells
celltype = c("Naïve CD4+ T-cell","TSCM","CD4+ TCM", 'aCD4+ TCM','aNaïve CD4+ T-cell',"IFN+ CD4+ T-cell",
             "cTreg","eTreg","Th2","TRM","Th17",'Th','CCL5+ CD4+ T-cell',
             'CD8+ TCM',"Naïve CD8+ T-cell","GZMK+ CD8+ T-cell","GZMK+ CCL3+ CCL4+ CD8+ T-cell",
             "GZMB+ GZMH+ CD8+ T-cell","TEMRA",'GZMK+ IL7R+ CD8+ T-cell')

celltype_factor = c("Naïve_CD4_T_cell","TSCM","CD4_TCM", 'aCD4_TCM','aNaïve_CD4_T_cell',"IFN_CD4_T_cell",
                    "cTreg","eTreg","Th2","TRM","Th17",'Th','CCL5+ CD4+ T-cell',
                    'CD8+ TCM',"Naïve_CD8_T_cell","GZMK_CD8_T_cell","GZMK_CCL3_CCL4_CD8_T_cell",
                    'GZMK_IL7R_CD8_T_cell',"TEMRA","GZMB_GZMH_CD8_T_cell")
celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]
unique(Idents(data_run_subset_label))[!(unique(Idents(data_run_subset_label)) %in% celltype)]

celltype_factor <- gsub("\\+","",celltype_factor)
celltype_factor<- gsub(" ","_",celltype_factor )
celltype_factor <- gsub("-","_",celltype_factor )

Ident_orig = unique(Idents(data_run_subset_label))
unique(Idents(data_run_subset_label))

data_input  =data_run_subset_label[ ,Idents(data_run_subset_label) %in% celltype]

Ident_list = as.character(Idents(data_input))
unique(Ident_list)

celltype[!(celltype %in% unique(Ident_list))]
Ident_orig[!(Ident_orig %in% celltype)]

unique(Ident_list)[!(unique(Ident_list) %in% celltype)]

#unique(Idents(data_input))[!( unique(Idents(data_input)) %in% celltype)]
mat = as.data.frame(t(data_input@assays[["RNA"]]@counts))


markers <- c("CD4","SELL","CCR7","JUN","FOSB","DUSP1","ISG15", "IFI6","IFI44L","IL7R","BCL2","CDK6","FOXP3","CTLA4","HLA-DRB1","CD74","GATA3","ITGB1","CRIP1","LGALS1","S100A10","S100A11","RORA","KLRB1","CCL5","GZMK","CD8A","CD8B","CCL3","CCL4","GZMH","GZMB","PRF1","GNLY","FCGR3A")
marker_counts <- mat[,colnames(mat) %in% markers]
rownames(marker_counts) <- gsub("\\.","-",rownames(marker_counts))
marker_counts$ident <- Ident_list
marker_counts$ident <- gsub("\\+","",marker_counts$ident )
marker_counts$ident <- gsub(" ","_",marker_counts$ident )
marker_counts$ident <- gsub("-","_",marker_counts$ident )

marker_counts$ident[!(marker_counts$ident %in% celltype_factor)]

marker_means <- aggregate(scale(marker_counts[, 1:(ncol(marker_counts)-1)]), list(marker_counts$ident), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means <- data.frame(t(marker_means))
marker_means$Marker <- rownames(marker_means)
marker_means_long <- gather(marker_means, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels=celltype_factor)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)
unique(marker_means_long$Cluster)

pdf(paste0("/disk2/Projects/EloRD/Paper/EloRD_Nivo_Oksana/", folder_name,"_T heatmap.pdf"),width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),
                       values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()
###############################################################
## Monocytes

labels = unique(Idents(data_run_subset_label))
celltype <-as.character( c("SELL+ CD14+ Mono","HLA-DR-low SELL+ CD14+ Mono","RNASE2-high SELL+ CD14+ Mono","TGFb1+ CD14+ Mono","CCL3+ CD14+ Mono",
                "IFN+ CD14+ Mono","CD14+ CD16+ Mono","CD16+ Mono","IFN+ CD16+ Mono",'cDC1','cDC2',
                "Neutrophil","MDPC"))

#labels = gsub("[\r\n]", "", labels)
#celltype = gsub("[\r\n]", "", celltype)



celltype[!(celltype %in% labels)]
labels[!(labels %in% celltype)]

data_input  =data_run_subset_label[ ,(Idents(data_run_subset_label)) %in% celltype]
unique(Idents(data_run_subset_label))
unique(Idents(data_input))

Ident_list = as.character(Idents(data_input))
Ident_list <- gsub("Mono","Monocytes",Ident_list )
celltype_factor <- gsub("Mono","Monocytes",celltype )


celltype[!(celltype %in% unique(Idents(data_run_subset_label)))]
mat = as.data.frame(t(data_input@assays[["RNA"]]@data))

markers <- c("CD14","SELL","S100A8", "S100A9", "S100A12",
             "HLA-DRA","HLA-DRB1","CD74","RNASE2","TIMP1","TGFB1","CD300E","LGALS3","LGALS2","CCL3","CCL3L3","IL1B","CXCL2","CXCL8",
             "MX1","ISG15","IFI6","IFI44L","FCGR3A","CSF1R","MS4A7","CLEC9A","CLEC10A","FCER1A",
             "CD1C","MPO","ELANE","PRTN3","AZU1","PCNA","MKI67","TOP2A")
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

pdf(paste0('/disk2/Projects/EloRD/Paper/EloRD_Nivo_Oksana/',folder_name,'_Mono heatmap.pdf'),width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + 
  scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp))))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()


## RSP Mono
umap <- read.csv("/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2/Data/data_Mono_DC.csv", 
                 row.names = 1)
data_input  =data_run_subset_label[ ,(Idents(data_run_subset_label)) %in% celltype]
mat_RSP = as.data.frame(t(data_input@assays[["RNA"]]@data))

mat_RSP_tmp <- read.delim("/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2/Data/matrix_Mono_DC.tsv")

mat_RSP_tmp <- read_tsv("/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2/Data/matrix_Mono_DC.tsv")
mat_RSP_tmp <- data.frame(t(mat_RSP_tmp))
rownames(mat_RSP_tmp) <- gsub("\\.","-",rownames(mat_RSP_tmp))

umap_plot <- merge(umap, mat_RSP, by = "row.names")

umap_plot$ident <- as.character(umap_plot$ident)
umap_plot$ident <- as.character(umap_plot$label)

umap_plot$label <- ifelse(umap_plot$ident == "4", "CD16+ Mono",
                          ifelse(umap_plot$ident == "28", "IFN+ CD16+ Mono",
                                 ifelse(umap_plot$ident %in% c("2","19"), "SELL+ CD14+ Mono",
                                        ifelse(umap_plot$ident == "11", "IFN+ CD14+ Mono",
                                               ifelse(umap_plot$ident == "12", "CCL3+ CD14+ Mono",
                                                      ifelse(umap_plot$ident == "0", "TGFb1+ CD14+ Mono",
                                                             ifelse(umap_plot$ident == "9", "HLA-DR-low SELL+ CD14+ Mono",
                                                                    ifelse(umap_plot$ident == "10", "RNASE2-high SELL+ CD14+ Mono",
                                                                           ifelse(umap_plot$ident == "6", "CD14+ CD16+ Mono",
                                                                                  ifelse(umap_plot$ident == "5", "cDC2",
                                                                                         ifelse(umap_plot$ident == "27", "cDC1",
                                                                                                ifelse(umap_plot$ident == "16", "Neutrophil",
                                                                                                       ifelse(umap_plot$ident == "8", "MDPC","Remove")))))))))))))

cell_types <- c("SELL+ CD14+ Mono","HLA-DR-low SELL+ CD14+ Mono","RNASE2-high SELL+ CD14+ Mono","TGFb1+ CD14+ Mono","CCL3+ CD14+ Mono","IFN+ CD14+ Mono","CD14+ CD16+ Mono","CD16+ Mono","IFN+ CD16+ Mono","cDC1","cDC2","Neutrophil","MDPC")

markers <- c("CD14","SELL","S100A8", "S100A9", "S100A12","HLA-DRA","HLA-DRB1","CD74","RNASE2","TIMP1","TGFB1","CD300E","LGALS3","LGALS2","CCL3","CCL3L3","IL1B","CXCL2","CXCL8","MX1","ISG15","IFI6","IFI44L","FCGR3A","CSF1R","MS4A7","CLEC9A","CLEC10A","FCER1A","CD1C","MPO","ELANE","PRTN3","AZU1","PCNA","MKI67","TOP2A")

colnames(umap_plot) <- gsub('\\.',"-",colnames(umap_plot))
marker_df <- umap_plot[!as.character(umap_plot$label) %in% "Remove",colnames(umap_plot) %in% c(markers,"label")]
marker_means <- aggregate(scale(marker_df[, 2:ncol(marker_df)]), list(marker_df$label), mean)
rownames(marker_means) <- marker_means$Group.1
marker_means$Group.1  <- NULL
marker_means_t <- data.frame(t(marker_means))
colnames(marker_means_t) <- rownames(marker_means)
marker_means_t$Marker <- rownames(marker_means_t)
marker_means_long <- gather(marker_means_t, "Cluster", "Mean_Exp",-Marker)
marker_means_long$Cluster <- factor(marker_means_long$Cluster, levels=cell_types)
marker_means_long$Marker <- factor(marker_means_long$Marker, levels=markers)

pdf(paste0('/disk2/Projects/EloRD/Paper/EloRD_Nivo_Oksana/',folder_name,'Full_Run_MonoDC_Cell_heatmap.pdf'),width=12)
ggplot(marker_means_long) +geom_tile(aes(Marker,Cluster,fill=Mean_Exp)) + scale_fill_gradientn(colors = c("lightblue","white","orangered3"),values=scales::rescale(c(min(marker_means_long$Mean_Exp),0,max(marker_means_long$Mean_Exp)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme(axis.text.x = element_text(angle = 65, hjust=1)) + xlab("") + labs(fill="Mean Z-score") +theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=14)) +ylab("")
dev.off()


## Comparing Sylvia and RSP
mat == mat_RSP_tmp

marker_means_long_RSP == marker_means_long

