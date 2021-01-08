#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 03:14:52 2020

@author: sujwary
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 15:36:59 2020

@author: sujwary
"""

import os
import pandas as pd
import os.path
from os import path

filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]
#metaData = metaData[metaData['Run']== 1]

base_bam = '/disk2/Projects/EloRD/Data/Bam/'
output = '/disk2/Projects/EloRD/Data/coverage/'
os.chdir(output)


sample_name = 'GL1502BM'


## Sort and reindex bams
for i in range(1,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'samtools sort ' +  base_bam + sample_name + '.bam ' 
    command = command + '> ' + base_bam + sample_name + '_sorted.bam ' 
    print(command)
    os.system(command)

    command = 'samtools index ' +  base_bam + sample_name + '_sorted.bam ' 
    print(command)
    os.system(command)
    
    
for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'bedtools coverage -a ' 
    command = command + base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf '
    command = command + '-b ' + base_bam + sample_name + '_sorted.bam '
    command = command + '> ' + sample_name + '_coverage.bed'
    
    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if (path.exists(base_bam + sample_name + '_sorted.bam') 
        and  path.exists(base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf') 
        and not (path.exists(output + sample_name + '_coverage.bed')) ):
        print(command)
        os.system(command)
        
        
## Old
for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    print(sample_name)
    
    command = 'bedtools coverage -a ' 
    command = command + base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf '
    command = command + '-b ' + base_bam + sample_name + '.bam '
    command = command + '> ' + sample_name + '_coverage.bed'


        
    #runCommand = path.exists(base + sample_name + '.vcf') and not (path.exists(base + sample_name + '.vep.maf'))
    if (path.exists(base_bam + sample_name + '.bam') 
        and  path.exists(base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf') 
        and not (path.exists(output + sample_name + '_coverage.bed')) ):
        print(command)
        os.system(command)


# Check which files are missing      

for i in range(0,(metaData.shape[0])):
    sample_name = metaData['Sample'].iloc[i]
    #print(sample_name)

    if not (path.exists(base_bam + sample_name + '.bam')):
        #print(sample_name)
        print('')
        
    if not path.exists(base_bam + sample_name + '_demux_data_harmony/' + 'souporcell_merged_sorted_vcf.vep.vcf'):
        print(sample_name)
        print('')
        
