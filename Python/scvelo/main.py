# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print("Hi, {0}".format(name))  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

import os
import anndata
from pylab import rcParams
import matplotlib as mpl
import matplotlib
from matplotlib import rcParams

import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

base_input = '/home/sujwary/Desktop/scRNA/Output/Harmony/AllSamples/Batch_Sample_Kit//Subcluster/NK/Cluster/PCA30/res3/Data/'

sample = anndata.read_loom(base_input + "data.loom")

sample_obs = pd.read_csv(base_input + "cellID_obs.csv")
umap_cord = pd.read_csv(base_input + "cell_embeddings.csv")
cell_clusters = pd.read_csv(base_input + "clusters.csv")
