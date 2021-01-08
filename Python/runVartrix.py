#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:18:07 2020

@author: sujwary
"""
import os
import pandas as pd
import os.path
from os import path

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]

base_bam = '/disk2/Projects/EloRD/Data/Bam/'


for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = './vartrix -v ' 
    command = command  + base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf '
    command = command + '-b ' + base_bam + sample_name + '.bam '
    command = command + '-f ' + '/disk2/Projects/EloRD/Data/ReferenceData/refdata-cellranger-GRCh38-1.2.0.fasta '
    command = command + '-c ' + base_bam + sample_name + '_out_cell_barcodes_filter.csv '
    command = command + '-s ' + 'coverage ' 
    command = command + '--ref-matrix' + base_bam + 'vartrix/' + sample_name + '_coverage.mtx '
    command = command + '-o ' + + base_bam + 'vartrix/' + sample_name + '_output.mtx '
    

    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if (path.exists(base_bam + sample_name + '.bam') 
        and  path.exists(base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf') 
        and not (path.exists(base_bam + 'vartrix/' + sample_name + '_coverage.mtx ' ) )):
        print(command)
        os.system(command)
        
  # Check which files are missing      
for i in range(0,(metaData.shape[0])):

    if (path.exists(base_bam + sample_name + '.bam') 
        and  path.exists(base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf') 
        and not (path.exists(base_bam + 'vartrix/' + sample_name + '_coverage.mtx ' ) )):
        print(sample_name)

