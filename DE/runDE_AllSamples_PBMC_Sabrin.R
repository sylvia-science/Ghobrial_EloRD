
data_harmony_run_label_subset$Sample_type[DE_input$'Sample Type' == "PBMC"] = 'PB'
data_harmony_run_label_subset$Sample_type[DE_input$'Sample Type' == "Bone Marrow"] = 'BM'


DE_input = data_harmony_run_label_subset

# Define the idents you want to compare
celltype = 'NK'
ident1 = paste0('baseline BM ',celltype)
ident2 = paste0('baseline PB ',celltype)


# Make DE_ident variable that has your desired idents
DE_input$DE_ident = paste0(DE_input$Treatment,' ', DE_input$Sample_type, ' ', Idents(DE_input))


# Subset to only those two idents
DE_input = DE_input[,DE_input$DE_ident %in% c(ident1,ident2)]
ncol(DE_input)
unique(DE_input$DE_ident)

IdentPerSample = table(DE_input$sample )

sample_list = names(IdentPerSample)
sample_keep = sample_list[IdentPerSample >= 20]

# Get variables that you will put into design matrix
kit = factor(DE_input$kit)
ident = factor(as.character(Idents(DE_input)))
DE_input$ident = ident
DE_ident = factor(DE_input$DE_ident)
Patient = factor(DE_input$`Patient Number`)
DE_input$Patient = Patient
Treatment = factor(DE_input$Treatment)


# Make formula and design matrix
formula = ~DE_ident + kit + Patient

design <- model.matrix(formula)
colnames(design) <- gsub("DE_ident", "", colnames(design))
print(colnames(design))

# If design matrix isn't full rank, remove columns that are redundant
print(is.fullrank(design))
design_fr = fullRank(design)
is.fullrank(design_fr)
colnames_old = colnames(design)
colnames_fr = colnames(design_fr)

# Make contrast matrix that says you want to compare ident2 to ident1
int2 = match(ident2,colnames(design_fr) )
con <- integer(ncol(design_fr))
con[int2] <- 1

# Find genes that are well expressed which you will keep
keep <- filterByExpr(DE_input@assays[["RNA"]]@counts, group=DE_input$DE_ident, 
                     min.count = 1,min.total.count=10, 
                     large.n = 10,min.prop = 0.1)

result_edgeR = runEdgeR(DE_input,design_fr, contrast = con, keep = keep,
                        folder_output = filepath_cluster,
                        subfolder = subfolder)
