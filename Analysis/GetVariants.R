
library(Rsamtools)
library(Matrix)
library(bedr)
library(vcfR)
library(maftools)
library(readxl)
library(numbers)

source('~/Desktop/scRNA/Code/Integrate All/LoadHarmonyData.R')


base = '/disk2/Projects/EloRD/Data/Bam/'
sample = 'GL1305BM'

# Seurat
filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)
#metaData = metaData[metaData$Run== 1 || metaData$`Sample Type` == 'PBMC',]
#metaData = metaData[metaData$Run== 1,]

filename_sampleParam <- paste0('/home/sujwary/Desktop/scRNA/Data/sample_','Combine','_parameters.xlsx')
sampleParam_combine <- read_excel(filename_sampleParam)

bed <- read.table("/disk2/Projects/EloRD/Data/coverage/GL1003BM_coverage.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

colnames(bed) = c('chrom','chromStart','chromEnd','ref','alt','score','strand',
                  'thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts')

max(bed)

data_label = data_harmony_run_label
sample_list = unique(data_label$sample )

#sample_list = metaData$Sample

plot = DimPlot(data_label,pt.size = 0.7, reduction = "umap",label = TRUE,
               label.size = 6)
print(plot)



maf_genes_unique = unique(maf_all_meta$Hugo_Symbol)
maf_genes_unique  [ !(maf_genes_unique %in% rownames(data_label) )]
gene_list = c('OR37G','OR9-8','PTK5','GTK')
gene_list [gene_list %in% rownames(data_label)]


tmp = maf_all_meta[maf_all_meta$Hugo_Symbol == 'FRK',]


tmp = maf_all_meta[maf_all_meta$Hugo_Symbol == 'OR13C3',]

########################
## Save cell barcodes
########################
for (i in 1:length(sample_list)){
  sample = sample_list[i]
  print(sample)
  #data_label = data[[1]]
  #filepath_cluster = data[[2]]
  
  #data_label = data_label[,data_label$sample == sample]
  #barcode_harmony = colnames(data_label)
  #barcode_harmony <-gsub("_.*","",barcode_harmony)
  
  filename = paste("/home/sujwary/Desktop/scRNA/Data/",sample,"_raw_feature_bc_matrix.h5",sep = "")
  #data_sample = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  #barcode = colnames(data_sample)
  #barcode = unique(c(barcode,barcode_harmony))
  print(filename)
  #file = paste0(base, sample,'_out_cell_barcodes.csv')
  #write.table(barcode,file, sep = ',', row.names = F, col.names = F, quote = F)

  data_label_subset = data_label[,data_label$sample == sample]
  barcode_filter = sub("_.*", "", colnames(data_label_subset))
  
  file = paste0(base, sample,'_out_cell_barcodes_filter.csv')
  write.table(barcode_filter,file, sep = ',', row.names = F, col.names = F, quote = F)
  
  
}
###################################
## Load files and save maf data
###################################
str = '_harmony'
str = ''
for (i in 1:length(sample_list)){
  #sample_name = sampleParam$Sample[i]
  sample = sample_list[i]
  print(sample)

  #data_label = data[[1]]
  #filepath_cluster = data[[2]]
  if ( !file.exists(paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv'))){
    next
  }

  # barcodes
 
  barcodes = read.csv(paste0(base, sample, '_out_cell_barcodes',str,'.csv'), header = F)
  barcodes = as.character(barcodes$V1)
  length(barcodes)
  
  #barcode_filter = read.csv(paste0(base,sample,"_out_cell_barcodes_cellranger.csv"), header = F)
  #barcode_filter = as.character(barcode_filter$V1)
  #barcodes = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
  data_label_subset = data_label[,data_label$sample == sample]
  barcode_filter = sub("_.*", "", colnames(data_label_subset))

  length(barcode_filter)
  ## MMTX
  #str = '_harmony'
  ref = readMM(paste0(base, sample,'_demux_data',str,'/ref.mtx'))
  alt = readMM(paste0(base, sample,'_demux_data',str,'/alt.mtx'))
  dim(ref)
  dim(alt)
  colnames(ref) = barcodes
  colnames(alt) = barcodes
  
  ref = ref[,colnames(ref)  %in% barcode_filter]
  alt = alt[,colnames(alt)  %in% barcode_filter]
  ## VCF

  file_vcf = paste0(base, sample,
                    '_demux_data', str,'/',
                    'souporcell_merged_sorted_vcf.vcf')
  souporcell_merged_sorted = read.vcf(file_vcf)
  
  nrow(souporcell_merged_sorted[["vcf"]])

  #nrow(cluster_genotypes[["vcf"]])
  #nrow(cluster_genotypes_vep[["vcf"]])
  
  vcf = souporcell_merged_sorted[["vcf"]]
  
  
  
  maf = read.maf(maf = paste0(base, sample,'_demux_data', str,'/', sample,'.vep.maf'),
                 useAll = TRUE)
  
  
  #cluster_list = levels(unique(Idents(data_label)))
  

  #vcf_subset = vcf[rowSums(alt)!= 0 || rowSums(ref)!= 0,]
  
  #maf_subset = maf@data[maf@data[["Start_Position"]] %in% vcf_subset$POS]
  #maf_silent_subset = maf@maf.silent[maf@maf.silent[["Start_Position"]] %in% vcf_subset$POS]
  
  colnames = c('Barcode','Sample','CellType','Alt','Ref',colnames(maf@data))
  maf_output = setNames(data.frame(matrix(ncol = length(colnames), nrow=0)), colnames)
  
  ident_list = Idents(data_label_subset)
  names(ident_list) <-gsub("_.*","",names(ident_list) )
  
  cnt = 1
  for(bc in colnames(ref)){
    
    if (mod(cnt,100) == 0){
      print(cnt)
      #print(nrow(maf_output))
    }
    cnt = cnt + 1
    cellType = as.character(ident_list[names(ident_list) == bc])
    #cellType = ''
    alt_list = alt[,bc]
    ref_list = ref[,bc]
    pos_list = vcf$POS
    
    maf_subset = maf@data
    maf_subset = maf_subset[maf_subset$Start_Position %in% pos_list]
    maf_pos = maf_subset$Start_Position
    pos_in_maf = pos_list %in% maf_pos
    alt_list = alt_list[pos_in_maf]
    ref_list = ref_list[pos_in_maf]
    pos_list = pos_list[pos_in_maf]
    
    match_pos_in_maf = match(pos_list,maf_pos)
    alt_list = alt_list[match_pos_in_maf]
    ref_list = ref_list[match_pos_in_maf]
    pos_list = pos_list[match_pos_in_maf]
    
    maf_subset$Barcode = bc
    maf_subset$Alt = alt_list
    maf_subset$Ref = ref_list
    maf_subset$CellType = cellType
    
    maf_subset = maf_subset[maf_subset$Alt>0 | maf_subset$Ref>0,]
    maf_subset$Sample = sample

    maf_output = rbind(maf_output,maf_subset)
      
  }
  
  mask = as.logical(colSums(is.na(maf_output)) != nrow(maf_output))
  col_list = colnames(maf_output)
  col_list = col_list[mask]
  maf_output = data.frame(maf_output)
  
  maf_output = maf_output[, col_list]
  
  maf_output = maf_output[maf_output$Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site'),]
  
  maf_output
  file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  write.csv(maf_output, file = file)

  
  #file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  #write.csv(output_maf, file = file)
}

########################
## Combine maf files
########################
sample = sample_list[1]
print(sample)
file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
maf_output_1 = read.csv( file = file)
mask = as.logical(colSums(is.na(maf_output_1)) != nrow(maf_output_1))
col_list = colnames(maf_output_1)
col_list = col_list[mask]
maf_output_1 = data.frame(maf_output_1)

maf_output_1 = maf_output_1[, col_list]


colnames = c('Barcode','Sample','CellType','Alt','Ref',colnames(maf_output_1),"CLIN_SIG", "PUBMED")
maf_all = setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)

for (i in 1:length(sample_list)){
  sample = sample_list[i]
  print(sample)
  file = paste0('/disk2/Projects/EloRD/Output/MafResults/',sample,'.csv')
  if (file.exists(file)){
    
    maf_output = read.csv( file = file)
    
    mask = as.logical(colSums(is.na(maf_output)) != nrow(maf_output))
    col_list = colnames(maf_output)
    col_list = col_list[mask]
    maf_output = data.frame(maf_output)
    

    Missing1 <- setdiff(colnames(maf_all), colnames(maf_output))  # Find names of missing columns
    maf_output[Missing1] <- NA
    
    Missing2 <- setdiff(colnames(maf_output), colnames(maf_all))  # Find names of missing columns
    maf_all[,Missing2] <- NA
    
    if (nrow(maf_output) > 0){
      maf_output$Sample = sample
      maf_all = rbind(maf_all,maf_output)
    }
  }
}

#####################
## Filter maf files
#####################

remove_cells = c('Remove','Mono/T-Cell DBL','Mono/CD8+ T Cell DBL')
maf_all = maf_all[!(maf_all$CellType %in% remove_cells),]
maf_all = maf_all[maf_all$Alt != 0,]
maf_all$VAF = maf_all$Alt/(maf_all$Alt + maf_all$Ref)



#patient_list = unique(metaData$`Patient Number`[metaData$`Sample Type` == 'PBMC'])

#metaData_subset = metaData[metaData$`Patient Number` %in% patient_list,]

#maf_all = maf_all[maf_all$VAF < 0.5,]
maf_all_meta = merge(maf_all, metaData, by="Sample")
maf_all_meta$CellType = as.character(maf_all_meta$CellType)
maf_all_meta$CellType[maf_all_meta$'Sample Type' == 'PBMC'] = 'PBMC'

sample_list = unique(maf_all_meta$Sample)
patient_list = unique(maf_all_meta$`Patient Number`)

pos_list = unique(maf_all_meta$Start_Position)
HGVSc_list = unique(maf_all_meta$HGVSc)

#tmp = apply(maf_all, 2, function(x) length(unique(x)))


maf_list <- vector("list", length = length(patient_list))
maf_list_keep <- vector("list", length = length(patient_list))

cnt = 1

colnames = c('Patient','CellType','Effect','PatientCellTypeEffect')
mutation_df = setNames(data.frame(matrix(ncol = length(colnames), nrow = 0)), colnames)

cellType_all_list = unique(maf_all_meta$CellType)
cellType_all_list = cellType_all_list[!(cellType_all_list %in% remove_cells)]

#maf_all_meta$keep = F

maf_matrix_all = (data.frame(matrix(ncol = length(cellType_all_list) + 1, 
                                nrow = 0)))
colnames(maf_matrix_all) = c(cellType_all_list, 'Patient Number')

for (patient in patient_list){
  print(patient)
  maf_patient = maf_all_meta[ maf_all_meta$`Patient Number` == patient,]
  cellType_list = unique(maf_patient$CellType)
  cellType_list = cellType_list[!(cellType_list %in% c('Remove'))]
  
  effect_list = unique(maf_patient$all_effects)
  maf_matrix = (data.frame(matrix(ncol = length(cellType_all_list), 
                                  nrow = length(effect_list))))
  colnames(maf_matrix) = cellType_all_list
  rownames(maf_matrix) = effect_list
  
  
  for (cellType in cellType_list){
    maf_sample_cellType = maf_patient[maf_patient$CellType == cellType,]
    
    effect_cellType = maf_sample_cellType$all_effects
    maf_matrix[,cellType] = effect_list %in% effect_cellType


    
    colnames = c('Patient','CellType','Effect','PatientCellTypeEffect')
    mutation_df_patient = setNames(data.frame(matrix(ncol = length(colnames), 
                                                     nrow = length(unique(maf_sample_cellType$all_effects)))), colnames)
    mutation_df_patient$Patient = patient
    mutation_df_patient$CellType = cellType
    mutation_df_patient$Effect = unique(maf_sample_cellType$all_effects)
    #mutation_df_patient$VAF = (maf_sample_cellType$VAF)

    mutation_df_patient$PatientCellTypeEffect = paste(mutation_df_patient$Patient,
                                                      mutation_df_patient$CellType, 
                                                      mutation_df_patient$Effect)
    
    mutation_df = rbind(mutation_df,mutation_df_patient)
  }
  
  effect_per_cell = rowSums(maf_matrix)
  effect_keep = effect_list[effect_per_cell< 4]
  maf_matrix_keep = maf_matrix[effect_per_cell< 4,]
  
  maf_list[[cnt]] = maf_matrix
  maf_list_keep[[cnt]] = maf_matrix_keep
  cnt = cnt + 1
  maf_matrix$`Patient Number` = patient
  maf_matrix_all = rbind(maf_matrix_all,maf_matrix)
  
  #file = paste0('/disk2/Projects/EloRD/Output/MafResults/Filter/Patient',patient,'.csv')
  #write.csv(maf_matrix_keep, file = file)
}

nrow(maf_matrix)
length(unique(mutation_df$PatientCellTypeEffect))

maf_all_meta_cell = maf_all_meta
mutation_df$keep = F
#maf_all_meta_cell$Effect = maf_all$all_effects

#maf_all_meta_cell[unique(maf_all_meta_cell$CellType)] <- NA

#col1 = colnames(maf_all_meta_cell)
#col2 = colnames(maf_matrix_all)

maf_matrix_all$all_effects = rownames(maf_matrix_all)
maf_matrix_all_subset = maf_matrix_all[,!(names(maf_matrix_all) %in% c('Patient Number','all_effects'))]
  
maf_matrix_all$keep = rowSums(maf_matrix_all_subset, na.rm = T)< 4
intersect(colnames(maf_all_meta_cell),colnames(maf_matrix_all))

maf_all_meta_cell = merge(maf_all_meta_cell,maf_matrix_all,  how='left')

intersect(colnames(maf_all_meta_cell),colnames(maf_matrix_all))


sum(maf_all_meta_cell$keep)
sum(!maf_all_meta_cell$keep)


keep = maf_all_meta_cell[maf_all_meta_cell$keep & 
                               maf_all_meta_cell$CellType != 'PBMC',]



pathName <-paste0('/disk2/Projects/EloRD/Output/MafResults/Plots/','allKeep.png')
png(file=pathName,width=600, height=600, res = 200)
title = ('All Keep')
plot = ggplot(keep, aes(x=VAF)) + geom_histogram() + 
  ggtitle(title)
print(plot)
dev.off()

not_keep = maf_all_meta_cell[!maf_all_meta_cell$keep & 
                               maf_all_meta_cell$CellType != 'PBMC',]
title = ('All Not Keep')
pathName <-paste0('/disk2/Projects/EloRD/Output/MafResults/Plots/','allNotKeep.png')
png(file=pathName,width=600, height=600, res = 200)
plot = ggplot(not_keep, aes(x=VAF)) + geom_histogram() + 
  ggtitle(title)
print(plot)
dev.off()


PBMC_Keep = maf_all_meta_cell[maf_all_meta_cell$keep & 
                                maf_all_meta_cell$PBMC == T & 
                                maf_all_meta_cell$CellType != 'PBMC',]
title = ('PBMC Keep')
pathName <-paste0('/disk2/Projects/EloRD/Output/MafResults/Plots/','PBMCKeep.png')
png(file=pathName,width=600, height=600, res = 200)
plot = ggplot(PBMC_Keep, aes(x=VAF)) + geom_histogram() + 
     ggtitle(title)
print(plot)
dev.off()
noPBMC_Keep = maf_all_meta_cell[maf_all_meta_cell$keep & 
                                   maf_all_meta_cell$PBMC == F & 
                                   maf_all_meta_cell$CellType != 'PBMC',]
title = ('No PBMC keep')
pathName <-paste0('/disk2/Projects/EloRD/Output/MafResults/Plots/','NoPBMCKeep.png')
png(file=pathName,width=600, height=600, res = 200)
plot = ggplot(noPBMC_Keep, aes(x=VAF)) + geom_histogram() + 
  ggtitle(title)
print(plot)
dev.off()

plot = ggplot(maf_all, aes(x=Alt)) + geom_histogram() + 
  ggtitle(title) + xlim(-1, 10)
print(plot)


maf_all

plot = ggplot(maf_all, aes(x=VAF)) + geom_histogram() + 
  ggtitle(title)
print(plot)

metaData$Patient = metaData$`Patient Number`
mutation_df = merge(mutation_df, metaData, by="Patient")

mutation_df = mutation_df[mutation_df$keep,]
mutation_df_SMM = mutation_df[mutation_df$`Diagnosis` == 'High Risk SMM',]
mutation_df_NBM = mutation_df[mutation_df$`Diagnosis` == 'NBM',]
  
variant_SMM = mutation_df_SMM[ !(mutation_df_SMM$Effect %in% mutation_df_NBM$Effect), ]


##########
## Analyze Bam files
##########

base = '/disk2/Projects/EloRD/Data/Bam/'
for (i in 1:length(sample_list)){
  sample = sample_list[i]
  print(sample)
  bamFile = paste0(base,sample,'.bam')
  sample_bam = readBAM(bamFile)
}