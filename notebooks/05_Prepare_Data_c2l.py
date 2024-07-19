#!/bin/bash
#######################################################################################################################
# Spatial transcriptomics unveils the in situ cellular and molecular hallmarks of the lung in fatal COVID-19
# Author: Carlos A. Garcia-Prieto
# Description: this script was used for the preparation of the processed Visium ST and HLCA core objects for spatial 
# mapping and deconvolution with cell2location (c2l).
# Usage: change data directory (dir) variable and execute the script step-by-step
# We followed cell2location and single-cell best practices spatial deconvultion tutorials
# For more information please see https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
# and https://www.sc-best-practices.org/spatial/deconvolution.html#cell2location
#######################################################################################################################

###Import modules
import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd
import pathlib
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

pd.set_option('display.max_columns', 50)
sc.logging.print_header()
print(f"squidpy=={sq.__version__}")

####Read processed Visium ST and HLCA data
adata_vis = ad.read("/mnt/beegfs/cgarcia/Spatial/COVID19/cell2location/HLCA_publication/HLCA/Concatenation_adata_hdg_6000.h5ad") #Concatenated VisiumST filtered data
adata_ref_core = ad.read("/mnt/beegfs/cgarcia/Spatial/COVID19/cell2location/HLCA_publication/HLCA/HLCA_core_parenchyma_filter.h5ad") #Filtered HLCA core object with parenchyma cells

###Set anndata reference raw counts as adata_ref.X
adata_ref_core.layers['normalized_counts'] = adata_ref_core.X.copy()
adata_ref_core.X = adata_ref_core.raw.X.copy()
adata_ref_core.layers['raw_counts'] = adata_ref_core.X.copy()

###Set data directory
dir = "/mnt/beegfs/cgarcia/Spatial/COVID19/cell2location/HLCA_publication/HLCA/"

###UMAP with cell types
with plt.rc_context():  # Use this to set figure params like size and dpi
    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.umap(adata_ref_core, color=["ann_finest_level"], wspace=0.4, show=False, use_raw=False)
    plt.savefig(f"{dir}adata_ref_clusters.png",dpi=300, format="png")
    plt.close()

###Remove mitochondrial genes
# find mitochondrial (MT) genes
adata_ref_core.var["MT_gene"] = [
    gene.startswith("MT-") for gene in adata_ref_core.var["feature_name"]
]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_ref_core.obsm["MT"] = adata_ref_core[:, adata_ref_core.var["MT_gene"].values].X.toarray()
adata_ref_core = adata_ref_core[:, ~adata_ref_core.var["MT_gene"].values]

###Subset both datasets to the same gene set which is the baseline for the mapping between the single cell and spatial data.
shared_features = [
    feature for feature in adata_vis.var_names if feature in adata_ref_core.var_names
]
adata_ref_core = adata_ref_core[:, shared_features].copy()
adata_vis = adata_vis[:, shared_features].copy()

###Fitting the reference model
#Gene selection
import cell2location as c2l

selected = c2l.utils.filtering.filter_genes(
    adata_ref_core, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
)

###Filter genes
adata_ref_core = adata_ref_core[:, selected].copy()
adata_vis = adata_vis[:, selected].copy()

###Save c2l preprocessed objects
adata_ref_core.write_h5ad(f"{dir}adata_ref_core_c2l_preproc.h5ad", compression="gzip")
adata_vis.write_h5ad(f"{dir}adata_vis_c2l_preproc.h5ad", compression="gzip")
