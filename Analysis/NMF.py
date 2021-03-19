import signatureanalyzer as sa
import pandas as pd

# ---------------------
# RUN SIGNATURE ANALYZER
# ---------------------

#celltype = 'NK_RemoveRiboFeatures'
#celltype = 'TCell_NK_Mono'
#celltype = ''
#path = "/home/sujwary/Desktop/scRNA/Data/NMF/Harmony_AllSamples_Sample_Kit_" + celltype
base = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/T Cell/Cluster/PCA30/res3.5/Data/'
path = base + 'matrix_T Cell_varFeature.tsv'

base = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/NK/Cluster/PCA30/res2/Data/'
path = base + 'matrix_NK_varFeature.tsv'

base = '/disk2/Projects/EloRD_Nivo_Oksana_PBMC_noNPBMC/Harmony/Batch_Sample_kit/Subcluster/Mono_DC/Cluster/PCA30/res2/Data/'
path = base + 'matrix_Mono_DC_varFeature.tsv'
input_matrix = pd.read_csv(path + "", sep='\t')

sa.run_matrix(matrix=input_matrix, outdir= base + 'phi1_alpha10/', nruns=100,
              verbose=True,plot_results=False,K0=40,tolerance = 0.0002,
              objective='gaussian', max_iter=30000, phi=1.0, a=10.0, 
              prior_on_W = 'L2',prior_on_H = 'L2',cuda_int = 0)
# tolerance=2e-04 
# k = 20 ceiling

# Priors L2




text_file = open("/home/sujwary/Desktop/scRNA/Data/NMF/PrePostEOTNBM_MT15_All_phi1_alpha10/Output.txt", "w")
text_file.write("Finished")
text_file.close()



# ---------------------
# LOADING RESULTS
# ---------------------
import pandas as pd

H = pd.read_hdf('nmf_output.h5', 'H')
W = pd.read_hdf('nmf_output.h5', 'W')
Hraw = pd.read_hdf('nmf_output.h5', 'Hraw')
Wraw = pd.read_hdf('nmf_output.h5', 'Wraw')
feature_signatures = pd.read_hdf('nmf_output.h5', 'signatures')
markers = pd.read_hdf('nmf_output.h5', 'markers')
cosine = pd.read_hdf('nmf_output.h5', 'cosine')
log = pd.read_hdf('nmf_output.h5', 'log')

# Output for each run may be found at...
Hrun1 = pd.read_hdf('nmf_output.h5', 'run1/H')
Wrun1 = pd.read_hdf('nmf_output.h5', 'run1/W')
# etc...

# Aggregate output information for each run
aggr = pd.read_hdf('nmf_output.h5', 'aggr')

# ---------------------
# PLOTTING
# ---------------------
sa.pl.marker_heatmap(...)
sa.pl.signature_barplot(...)
sa.pl.stacked_bar(...)
sa.pl.k_dist(...)
sa.pl.consensus_matrix(...)