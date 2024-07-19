#!/bin/bash
#######################################################################################################################
# Spatial transcriptomics unveils the in situ cellular and molecular hallmarks of the lung in fatal COVID-19
# Author: Carlos A. Garcia-Prieto
# Description: this script was used to fit the HLCA reference model which estimates the reference cell type signature.
# Cell2location estimates the signatures from the reference data using a Negative Binomial regression model which can 
# also account for batch effects and additional categorical covariates. 
# Usage: change results_folder variable and execute the script in a GPU cluster.
# We followed cell2location and single-cell best practices spatial deconvultion tutorials
# For more information please see https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
# and https://www.sc-best-practices.org/spatial/deconvolution.html#cell2location
#######################################################################################################################

###Import modules
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
#Set results folders
results_folder = "/mnt/beegfs/cgarcia/Spatial/COVID19/cell2location/HLCA_publication/HLCA/"

ref_run_name = f'{results_folder}reference_signatures_finest'
run_name = f'{results_folder}cell2location_map_finest'

########################################################################
#Estimation of reference cell type signatures (NB regression) 
##The signatures are estimated from scRNA-seq data, accounting for batch effect, using a Negative binomial regression model.

adata_ref = ad.read(f"{results_folder}adata_ref_core_c2l_preproc.h5ad")
adata_ref

##Delete unnecessary raw slot
del adata_ref.raw

#Prepare anndata for the regression model
#We are passing sample variable as batch and assay, donor-id, tissue_sampling_method and tissue_dissociation_protocol variables as additional covariates to account for batch effects.
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, 
                        # 10X reaction / sample / batch
                        batch_key='sample', 
                        # cell type, covariate used for constructing signatures
                        labels_key='ann_finest_level', 
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=["assay", "donor_id", "tissue_sampling_method", "tissue_dissociation_protocol"]
                       )
#Create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref) 

#View anndata_setup as a sanity check
mod.view_anndata_setup()

########################################################################
#Training model.
##Now we train the model to estimate the reference cell type signatures.
##Note that to achieve convergence on your data (=to get stabilization of the loss) you may need to increase max_epochs=250.
##Also note that here we are using batch_size=2500 which is much larger than scvi-tools default and perform training on all cells in the data (train_size=1) - both parameters are defaults.

mod.train(max_epochs=250, use_gpu=True)

#In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

#Save model
mod.save(f"{ref_run_name}", overwrite=True)

#Save anndata object with results
adata_file = f"{ref_run_name}/adata_ref_finest_post_prob.h5ad"
adata_ref.write(adata_file)

#Assess model training
with mpl.rc_context({'axes.facecolor':  'white','figure.figsize': [6, 4]}):
    mod.plot_history(20)
    plt.savefig(f"{ref_run_name}/NBregression_model_training.png",dpi=300, format="png",pad_inches=0.2,bbox_inches="tight")
    plt.close()

#with mpl.rc_context({'axes.facecolor':  'white','figure.figsize': [6, 6]}):
#    mod.plot_QC()
#    plt.savefig(f"{ref_run_name}/NBregression_model_training_QC.png",dpi=300, format="png",pad_inches=0.2,bbox_inches="tight")
#    plt.close()
