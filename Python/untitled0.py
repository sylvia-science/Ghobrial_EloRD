#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 22:33:54 2020

@author: sujwary
"""
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import os

os.chdir('/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/NK_Remove_Nfeature2000/Cluster/PCA30/res1.2/Data/')


sample_one = anndata.read_loom("/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit/Subcluster/NK_Remove_Nfeature2000/Cluster/PCA30/res1.2/Data/data.loom")
 
sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters.csv")

sample_obs["x"]
#sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["x"])]


sample_one.obs.index
sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'Cell ID'})

umap = umap_cord

umap = umap.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = sample_one_index.merge(umap, on = "Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values


scv.pp.filter_and_normalize(sample_one)
scv.pp.moments(sample_one)
scv.tl.velocity(sample_one, mode = "stochastic")
scv.tl.velocity_graph(sample_one)
scv.pl.velocity_embedding(sample_one, basis='umap')


color = sample_one.uns['Cluster_colors']
