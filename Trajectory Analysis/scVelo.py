#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:01:44 2020

@author: sujwary
"""
import os
import anndata
import pandas 
import matplotlib as mpl
from matplotlib import rcParams


import scvelo as scv
import pandas as pd
import numpy as np

import re
cell_clusters_ordered_all = []

import seaborn as sns
from matplotlib.colors import rgb2hex
import matplotlib

cellType = 'NK_RemoveRiboFeatures'
filename_metaData = '/home/sujwary/Desktop/scRNA/Data/EloRD Meta.xlsx'
metaData = pd.read_excel(filename_metaData)
metaData = metaData[metaData['Run']== 1]
sample_name = metaData['Sample'].iloc[2]

#base_input = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/NK_Remove_Nfeature2000/Cluster/PCA20/res1.2/Data/'

if cellType == 'T-cell':
    base_input= '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/T Cell/Cluster/PCA30/res3.5/Data/'
    xlim = [-6,6]
    ylim = [-10,14]

if cellType == 'NK_Remove_Nfeature2000':
    base_input= '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/T Cell/Cluster/PCA30/res3.5/Data/'
    xlim = [-6,6]
    ylim = [-10,14]
    
if cellType == 'NK':
    base_input=     '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/NK/Cluster/PCA30/res3/Data/'
    xlim = [-10,12]
    ylim = [-6,6]    
    remove = [0, 8, 11, 17, 18, 19, 21, 22]
  
if cellType == 'NK_RemoveRiboFeatures':
    base_input= '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK_RemoveRiboFeatures/Cluster/PCA30/res3/Data/'
    xlim = [-6,10]
    ylim = [-10,10]
      
  
input_loom = '/disk2/Projects/EloRD/Data/velocyto_harmony_filter/'

file_list = os.listdir(input_loom) 

cell_clusters = pd.read_csv(base_input + "clusters.csv")
all_cluster= set(cell_clusters['x'] )

num_cluster = len(set(cell_clusters['x'] ))
palette = sns.color_palette(None, num_cluster)
palette = np.array([mpl.colors.to_hex(x) for x in palette])

# sample_filter = anndata.read_loom('/disk2/Projects/EloRD/Data/velocyto_cell_barcodes_filter/' + sample_name +  ".loom")


#for i in range(0,(metaData.shape[0])):
for i in range(0,(metaData.shape[0])):
#for i in range(0,3):
    print(i)
    sample_name = metaData['Sample'].iloc[i]

    print(sample_name)
    sample_obs = pd.read_csv(base_input + "cellID_obs.csv")
    umap_cord = pd.read_csv(base_input + "cell_embeddings.csv")
    cell_clusters = pd.read_csv(base_input + "clusters.csv")


    loom_file = matching = [s for s in file_list if sample_name in s]
    if(len(loom_file) == 0):
        continue
    sample = anndata.read_loom(input_loom + loom_file[0])
    sample.var_names_make_unique()
    
    sample_obs["x"] = [re.sub("-.*", '', i) for i in sample_obs["x"]]
    sample_obs['sample'] = [re.sub("_.*", '', i) for i in sample_obs["x"]]
    
    umap_cord['Cell ID'] = sample_obs["x"]
    umap_cord['sample'] = sample_obs["sample"]
    cell_clusters['Cell ID'] = sample_obs["x"]
    cell_clusters['sample'] = sample_obs["sample"]
    
    
    sample_obs = sample_obs[sample_obs['sample'] == sample_name]
    umap_cord = umap_cord[umap_cord['sample'] == sample_name]
    cell_clusters = cell_clusters[cell_clusters['sample'] == sample_name]

    if len(sample_obs) == 0:
        continue
    print(sample_obs)
    
    
    barcode_loom = sample.obs.index.tolist()
    barcode_loom = [str(x) for x in barcode_loom]
    barcode_loom = [re.sub('^.*?:', '', i) for i in barcode_loom]
    barcode_loom = [sample_name + '_' + x[:-1] for x in barcode_loom]
    sample.obs.index = barcode_loom
    print(barcode_loom[0])
    print(sample_obs["x"].iloc[0])

    len(barcode_loom)
    len(sample_obs["x"])


    print(sum(np.isin(barcode_loom,sample_obs["x"])))
    sample_subset = sample[np.isin(barcode_loom,sample_obs["x"])]
    # sample_filter_subset = sample_filter[np.isin(barcode_loom_filter,sample_obs["x"])]
    
    len(sample_subset.obs.index.tolist())
    # len(sample_filter_subset.obs.index.tolist())
    

    sample_index = pd.DataFrame(sample_subset.obs.index.tolist())
    sample_index = sample_index.rename(columns = {0:'Cell ID'})
    len(sample_index)
    # sample_index_filter = pd.DataFrame(sample_filter_subset.obs.index.tolist())
    # sample_index_filter = sample_index_filter.rename(columns = {0:'Cell ID'})
    
    
    #tmp1 = umap_cord['Cell ID'][1]
    #tmp2 = sample_index['Cell ID'][2]
    #tmp = np.isin(umap_cord['Cell ID'],sample_index['Cell ID'])
    
    umap_ordered = sample_index.merge(umap_cord, on = "Cell ID")
    cell_clusters_ordered = sample_index.merge(cell_clusters, on = "Cell ID")

    print(all(sample_subset.obs.index.tolist() == umap_ordered['Cell ID']))

    umap_ordered = umap_ordered.iloc[:,2:4]
    cell_clusters_ordered = cell_clusters_ordered.iloc[:,2]
    len(umap_cord)
    len(umap_ordered)
    
    # umap_ordered_filter = sample_index_filter.merge(umap_cord, on = "Cell ID")
    # umap_ordered_filter = umap_ordered_filter.iloc[:,1:]
    
    sample_subset.obsm['X_umap'] = umap_ordered.values

    
    sample_subset.obs['clusters'] = cell_clusters_ordered.values
    cell_clusters_ordered = cell_clusters_ordered[~np.isin(cell_clusters_ordered.values,remove)]
    sample_subset = sample_subset[ ~np.isin(sample_subset.obs['clusters'].values,remove)]

    

    
    cluster_val = (list(set(cell_clusters_ordered)))
    
    palette_subset = palette[[i in cluster_val for i in   all_cluster]]
    sample_subset.uns['clusters_colors'] =palette_subset

    #sample_all.uns['clusters_colors'] = np.array(palette)


    cell_clusters_ordered_all.extend (cell_clusters_ordered.values)

    
#=============================================================================
    if (i == 0):
        
        print('') 
        sample_all = sample_subset
        continue
    else:
        sample_all = sample_all.concatenate(sample_subset)
        continue
#=============================================================================
    if (len(sample_subset) < 31):
        print('')
        continue
    
    
    # sample_filter_subset.obsm['X_umap'] = umap_ordered_filter.values
    #sample_one.uns['Cluster_colors']
    print(len( umap_ordered.values))
    # print(len(umap_ordered_filter.values))
    
    os.chdir('/disk2/Projects/EloRD/Output/scVelo/'+ cellType +'/')
    
    scv.pp.filter_and_normalize(sample_subset)
    scv.pp.moments(sample_subset)
    scv.tl.velocity(sample_subset, mode = "stochastic")
    scv.tl.velocity_graph(sample_subset)
    
    save_file = sample_name + '_Velocity_umap'+ '_color' + '.pdf'
    scv.pl.velocity_embedding(sample_subset, basis='umap',
                              xlim = xlim,ylim = ylim,
                              fontsize  = 12,save= save_file)

    try:
        save_file = sample_name + '_Velocity_stream_umap'+ '_color' + '.png'
        scv.pl.velocity_embedding_stream(sample_subset, basis='umap',
                                         xlim = xlim,ylim = ylim,
                                         figsize = [14,10],save= save_file)
    
    except ValueError:
        print('error')
        continue
            
################################################


cell_clusters_ordered_all = np.array(cell_clusters_ordered_all)
sample_all.obs['clusters'] = cell_clusters_ordered_all


num_cluster = len(set(cell_clusters_ordered_all ))



palette = sns.color_palette(None, num_cluster)

palette = [matplotlib.colors.to_hex(x) for x in palette]

sample_all.uns['clusters_colors'] = np.array(palette)


os.chdir('/disk2/Projects/EloRD/Output/scVelo/'+ cellType +'/')
    
scv.pp.filter_and_normalize(sample_all)
scv.pp.moments(sample_all)
scv.tl.velocity(sample_all, mode = "stochastic")
scv.tl.velocity_graph(sample_all)
    
save_file = 'All' + '_Velocity_umap.pdf'
scv.pl.velocity_embedding(sample_all, basis='umap',
                          xlim = xlim,ylim = ylim,save= save_file)
save_file = 'All' + '_Velocity_stream_umap.png'
scv.pl.velocity_embedding_stream(sample_all, basis='umap',color='clusters',
                                 density = 0.3, min_mass = 3,add_polyfit = True,
                                 xlim = xlim,ylim = ylim,
                                 figsize = [14,10],save= save_file)


##

