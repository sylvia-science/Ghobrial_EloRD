#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 01:30:10 2020

@author: sujwary
"""
import os
import pandas as pd
import os.path
from os import path
import csv

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]
#metaData = metaData[metaData['Run']== 1]

base_bam = '/disk2/Projects/EloRD/Data/Bam/'
output = '/disk2/Projects/EloRD/Data/velocyto_harmony_filter/'
#os.chdir(output)


sample_name = 'GL1003BM'

for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    kit = metaData['10X kit'].iloc[i]
    print(sample_name)

    #command = 'samtools sort ' + base_bam + sample_name + '.bam '
    #command = command + '-o '  + base_bam + sample_name + '.sorted.bam '

    #if (path.exists(base_bam + sample_name + '.bam')):
    #    print(command)
    #    os.system(command)
    input_csv = base_bam + sample_name + '_out_cell_barcodes_filter.csv'
    output_tsv = base_bam + sample_name + '_out_cell_barcodes_filter.tsv'
    
    #if (path.exists(input_csv)):
    #    csv.writer(open(output_tsv, 'w+'), delimiter='\t').writerows(csv.reader(open(input_csv))) 
        
    if (kit == 'v3_1'):
        barcode_file = '/disk2/Projects/10XBarcodeWhitelist/3M-february-2018_V3.tsv'
        
    if (kit == 'v2'):
        barcode_file = '/disk2/Projects/10XBarcodeWhitelist/737K-august-2016_V2.tsv'
    
    barcode_file = '/disk2/Projects/EloRD/Data/Harmony_Barcodes/' + sample_name + '_barcode.tsv'
    
    command = 'velocyto run -b ' 
    command = command + barcode_file + ' '
    command = command + '-o ' + output + ' '
    command = command + '-m /disk2/Projects/EloRD/Data/ReferenceData/repeat_mask.gtf '
    command = command + base_bam + sample_name + '.bam '
    command = command + '/disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf'

   #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if (path.exists(base_bam + sample_name + '.bam')):
        print(command)
        os.system(command)
