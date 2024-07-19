[![DOI](https://badgen.net/static/bioRxiv/10.1101%2023.12.21.572330/red?.svg)](https://doi.org/10.1101/2024.07.03.601404)

# *Spatial transcriptomics unveils the in situ cellular and molecular hallmarks of the lung in fatal COVID-19*

Carlos A. Garcia-Prieto, Eva Musulen, Daniela Grases, Veronica Davalos, Gerardo Ferrer, Belén Pérez-Miés, Tamara Caniego-Casas, José Palacios, Xavier Saenz-Sardà, Elisabet Englund, Eduard Porta, Manel Esteller

*bioRxiv* DOI:[https://doi.org/10.1101/2024.07.03.601404]


## Description

This repository contains the Python code used to produce all the analyses and figures presented in the article *Spatial transcriptomics unveils the in situ cellular and molecular hallmarks of the lung in fatal COVID-19* (Garcia-Prieto, CA et al.).

The data discussed in this publication, including all Visium datasets, sample metadata, raw and processed data (Space Ranger output), have been deposited in NCBI’s Gene Expression Omnibus (GEO) and are accessible through GEO Series accession number [GSE271370](https://www.ncbi.nlm.nih.gov./geo/query/acc.cgi?acc=GSE271370). 

Visium Spatial Gene Expression libraries for Formalin Fixed Paraffin Embedded (FFPE) tissue samples were analyzed with spaceranger (v. 2.0.0) count pipeline from 10x Genomics. First, a manual fiducial alignment and tissue boundary identification, including manual selection of spots covering tissue regions, were performed for each single library FFPE sample using Loupe Browser (v. 6.0). A probe set reference file compatible with FFPE workflow and human reference genome GRCh38 were downloaded from 10x Genomics and used to map Visium gene expression libraries. 

We leveraged the Human Lung Cell Atlas [(HLCA)](https://doi.org/10.1038/s41591-023-02327-2) single-cell RNA sequencing (scRNA-seq) reference dataset and consensus cell type annotations for spatial mapping and annotation with [cell2location](https://cell2location.readthedocs.io/en/latest/)(v. 0.1.3).

The overview of the study design, including sample processing workflow and spatial transcriptomics data analysis pipeline, have been illustrated here: 


<p align="center"><img src="/doc/Graphical_abstract.png" height="440" width=800"></p>


## Content

* `notebooks/`: All Python notebooks and scripts used to produce all the analyses and figures   
* `doc/`: Contains schematic workflow overview

## Contact

For questions related to this repo and its content, please contact Carlos A. Garcia-Prieto (cgarcia@carrerasresearch.org)
