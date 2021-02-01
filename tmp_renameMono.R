mono <- read.csv("/home/sujwary/Desktop/Mono_SampleByCellType.csv", row.names = 1,check.names=FALSE)
#mono <- read.csv("/home/sujwary/Desktop/Mono_SampleByCellType.csv")
mono = colSums(mono)

mono_new = read.csv('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/Mono_DC/Cluster/PCA30/res1.6/Stats/Mono_SampleByCellType.csv',
  row.names = 1, check.names = F)

mono_new = colSums(mono_new)
old_colnames = (colnames(mono))
old_new_colnames = (colnames(mono_new))
RSP_colnames <- (c("CD14+ Mono/T-cell DBL","cDC2","CD16+ Mono","TGFb1+ CD14+ Mono", "SELL+ CD14+ Mono", "sDC", "sMono","prDC","CD14+ CD16+ Mono","dDC","CD14+ Mono/CD8+ T-cell DBL","dMono","dMIP1a+ Mono","cDC1", "DC/T-cell DBL","IFN+ Mono", "GMPC","MK","sCD14+ Mono","MIP1a+ CD14+ Mono", "Erythrocytes"))

mono_sort = sort((mono))
mono_new_sort = sort((mono_new))
mono_sort
mono_new_sort


RSP_colnames[names(mono) == 'CD14+ Mono' ]
mono[R_colnames == 'dMIP1a+ Mono' ]

mono_sort [ 8]

mono_new = mono_new[names(mono_new)!= '14']
print(sort(names(mono_new)))
print(sort(R_colnames))

for (i in 1:length(new_colnames)){
  print('i: ')
  print(i)
  print('Old labels')
  print(mono[i])
  print('New labels')
  print(mono_new[i])
  print('RSP labels')
  print(R_colnames[i])
  print('')
  
}


