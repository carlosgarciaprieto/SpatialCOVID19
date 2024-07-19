#!/bin/bash
#######################################################################################################################
# Spatial transcriptomics unveils the in situ cellular and molecular hallmarks of the lung in fatal COVID-19
# Author: Carlos A. Garcia-Prieto
# Description: this script was used to prepare and train cell2location model to perform the spatial mapping and to 
# estimate cell-type specific expression of every gene in the spatial data (needed for NCEM).
# Cell2location is a Bayesian model that estimates the absolute abundance of cell types at each location by decomposing 
# the spatial expression count matrix into a predefined set of reference cell type signatures. 
# Importantly, cell2location allows the joint modeling of multiple spatial datasets  (also known as batch integration). 
# Usage: change results_folder variable and execute the script in a GPU cluster.
# We followed cell2location and single-cell best practices spatial deconvultion tutorials
# For more information please see https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
# and https://www.sc-best-practices.org/spatial/deconvolution.html#cell2location
#######################################################################################################################

#Import modules
import sys
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib as mpl
import scvi
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True' # this line should go before importing cell2location
TF_CPP_MIN_LOG_LEVEL=0
import cell2location
import torch
torch.set_float32_matmul_precision('medium')
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs
pd.set_option('display.max_columns', 100)

########################################################################
#Set folders
results_folder = "/mnt/beegfs/cgarcia/Spatial/COVID19/cell2location/HLCA_publication/HLCA/"

ref_run_name = f'{results_folder}reference_signatures_finest'
run_name = f'{results_folder}cell2location_map_finest'

########################################################################
#Read anndata_ref model with reference cell type signatures
adata_file = f"{ref_run_name}/adata_ref_finest_post_prob.h5ad"
adata_ref = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

inf_aver.to_csv(f"{ref_run_name}/inf_aver.csv")

########################################################################
#Cell2location: spatial mapping 
##Find shared genes and prepare anndata. Subset both anndata and reference signatures:
##Prepare Cell2location training model
adata_vis = sc.read_h5ad(f"{results_folder}adata_vis_c2l_preproc.h5ad")

#Find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

intersect.shape

#Prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, layer="raw_counts", batch_key="sample")

#Create and train the model
#Hyperparameters were set as in Madissoon, E., Nat Genet 2023 paper 
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver, 
    # the expected average cell abundance: tissue-dependent 
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=20,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
) 
mod.view_anndata_setup()

mod.train(max_epochs=40000, 
          # train using full data (batch_size=None)
          batch_size=None, 
          # use all data points in training because 
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/ad_vis_post_distrib_finest.h5ad"
adata_vis.write(adata_file)

# Assess model training
with mpl.rc_context({'axes.facecolor':  'white','figure.figsize': [6, 4]}):
    mod.plot_history(1000)
    #plt.legend(labels=['full data training'])
    plt.savefig(f"{run_name}/cell2location_model_training.png",dpi=300, format="png",pad_inches=0.2,bbox_inches="tight")
    plt.close()

with mpl.rc_context({'axes.facecolor':  'white','figure.figsize': [6, 6]}):
    mod.plot_QC()
    plt.savefig(f"{run_name}/cell2location_model_training_QC.png",dpi=300, format="png",pad_inches=0.2,bbox_inches="tight")
    plt.close()

# Plot QC metrics
with mpl.rc_context({'axes.facecolor':  'white','figure.figsize': [4.5, 5]}):
    fig = mod.plot_spatial_QC_across_batches()
    plt.savefig(f"{run_name}/spatial_QC_batches.png",dpi=300, format="png",pad_inches=0.2,bbox_inches="tight")
    plt.close()

### Estimate cell-type specific expression of every gene in the spatial data (needed for NCEM)

# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples[f"post_sample_q05"], mod.adata_manager #mod.samples
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
adata_file = f"{run_name}/ad_vis_post_distrib_finest_cell_type_gene_expr.h5ad"
adata_vis.write(adata_file)

