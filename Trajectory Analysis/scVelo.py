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

base_input = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK/Cluster/PCA30/res3/Data/'

sample = anndata.read_loom(base_input + "data.loom")
 

sample_obs = pd.read_csv(base_input + "cellID_obs.csv")
umap_cord = pd.read_csv(base_input + "cell_embeddings.csv")
cell_clusters = pd.read_csv(base_input + "clusters.csv")



sample_index = pd.DataFrame(sample.obs.index)
sample_index = sample_index.rename(columns = {0:'Cell ID'})


umap_cord.rename(columns={"Unnamed: 0":"Cell ID"}, inplace=True)

umap_ordered = sample_index.merge(umap_cord,on="Cell ID")

umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values
#sample_one.uns['Cluster_colors']



scv.pp.filter_and_normalize(sample)
scv.pp.moments(sample)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)
scv.pl.velocity_embedding(sample, basis='umap')