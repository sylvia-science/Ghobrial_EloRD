library(readxl)
library(plyr)
library(xlsx)


filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = read_excel(filename_metaData)


filename_patient = '/home/sujwary/Desktop/scRNA/Data/Elo_Patient info_updated_010120.xlsx'
patientData = read_excel(filename_patient)
patientData = patientData[,c('Processed Date','Sample','Patient Initials','Patient #','FISH','Subtype','Gender','DOB')]



metaData_merge = (join(metaData,patientData,by = 'Sample'))
#metaData_merge[is.na(metaData_merge)] <- 'NA'
write.xlsx(metaData_merge, file = paste0('/home/sujwary/Desktop/scRNA/Data/','EloRD Meta_12-23-2020.xlsx'), row.names = FALSE)
