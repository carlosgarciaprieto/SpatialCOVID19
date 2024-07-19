#!/bin/bash
#Script to run spaceranger count pipeline using Space Ranger 2.0.0 for FFPE tissue samples
#Probe set reference file compatible with FFPE workflow and human reference genome GRCh38 were downloaded from 10x Genomics
#All Visium datasets, sample metadata, raw and processed data (including Space Ranger output), have been deposited in NCBI Gene Expression Omnibus (GEO) 
#and are accessible through GEO Series accession number GSE271370(https://www.ncbi.nlm.nih.gov./geo/query/acc.cgi?acc=GSE271370)

spaceranger count --id=L2_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L2/HJ3F3CCX2 --sample=L2 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide1area1multifoco.tif --slide=V11B18-010 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-010.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L2.json
spaceranger count --id=L5_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L5/HJ3F3CCX2 --sample=L5 --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide1area2(1).tif' --slide=V11B18-010 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-010.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L5.json
spaceranger count --id=L14_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L14/HJ3F3CCX2 --sample=L14 --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide1area4(1).tif' --slide=V11B18-010 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-010.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L14.json
spaceranger count --id=L19_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L19/HJ3F3CCX2 --sample=L19 --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide2area1(1).tif' --slide=V11B18-011 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-011.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L19.json
spaceranger count --id=L24_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L24/HJ3F3CCX2 --sample=L24 --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide2area2(1).tif' --slide=V11B18-011 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-011.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L24.json
spaceranger count --id=L3C_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L3C/HJ3F3CCX2 --sample=L3C --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide2area3(1).tif' --slide=V11B18-011 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-011.gpr --area=C1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L3C.json
spaceranger count --id=L14C_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00174897/HN00174897_10X_RawData_Outs/L14C/HJ3F3CCX2 --sample=L14C --image='/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/Slide2area4(1).tif' --slide=V11B18-011 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/V11B18-011.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Results/L14C.json
spaceranger count --id=L2C_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00183526/HN00183526_10X_RawData_Outs/L2C/HJJY3CCX2 --sample=L2C --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L2C.tif --slide=V12U07-375 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/V12U07-375.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L2C.json
spaceranger count --id=L11P_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00183526/HN00183526_10X_RawData_Outs/L11P/HJJY3CCX2 --sample=L11P --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L11P.tif --slide=V12U07-376 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/V12U07-376.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L11P.json
spaceranger count --id=L12P_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00183526/HN00183526_10X_RawData_Outs/L12P/HJJY3CCX2 --sample=L12P --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L12P.tif --slide=V12U07-376 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/V12U07-376.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/SegundoBatch/L12P.json
spaceranger count --id=COVID_HRC2_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC2/HV3VMDSX5 --sample=COVID_HRC2 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC2.tif --slide=V12U20-281 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-281.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC2.json
spaceranger count --id=COVID_HRC4_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC4/HV3VMDSX5 --sample=COVID_HRC4 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC4.tif --slide=V12U20-281 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-281.gpr --area=C1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC4.json
spaceranger count --id=COVID_HRC5_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC5/HV3VMDSX5 --sample=COVID_HRC5 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC5.tif --slide=V12U20-281 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-281.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC5.json
spaceranger count --id=COVID_HRC6_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC6/HV3VMDSX5 --sample=COVID_HRC6 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC6.tif --slide=V12U20-282 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-282.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC6.json
spaceranger count --id=COVID_HRC8_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC8/HV3VMDSX5 --sample=COVID_HRC8 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC8.tif --slide=V12U20-282 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-282.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC8.json
spaceranger count --id=COVID_HRC10_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC10/HV3VMDSX5 --sample=COVID_HRC10 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC10.tif --slide=V12U20-282 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-282.gpr --area=C1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC10.json
spaceranger count --id=COVID_HRC11_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC11/HV3VMDSX5 --sample=COVID_HRC11 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC11.tif --slide=V12U20-282 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-282.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC11.json
spaceranger count --id=COVID_HRC_12_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC_12/HV3VMDSX5 --sample=COVID_HRC_12 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC12.tif --slide=V12U20-283 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-283.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC12.json
spaceranger count --id=COVID_HRC13_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC13/HV3VMDSX5 --sample=COVID_HRC13 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC13.tif --slide=V12U20-283 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-283.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC13.json
spaceranger count --id=COVID_HRC16_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC16/HV3VMDSX5 --sample=COVID_HRC16 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC16.tif --slide=V12U20-283 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-283.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC16.json
spaceranger count --id=COVID_HRC17_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC17/HV3VMDSX5 --sample=COVID_HRC17 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC17.tif --slide=V12U20-284 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-284.gpr --area=A1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC17.json
spaceranger count --id=COVID_HRC18_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/COVID_HRC18/HV3VMDSX5 --sample=COVID_HRC18 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC18.tif --slide=V12U20-284 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-284.gpr --area=B1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/COVID_HRC18.json
spaceranger count --id=CONTROL_2_outs --transcriptome=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/refdata-gex-GRCh38-2020-A --probe-set=/gpfs/scratch/bsc59/bsc59829/Spatial/Downloads/Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv --fastqs=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/HN00192875/HN00192875_10X_RawData_Outs/CONTROL_2/HV3VMDSX5 --sample=CONTROL_2 --image=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/CONTROL2.tif --slide=V12U20-284 --slidefile=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Slide_files/V12U20-284.gpr --area=D1 --loupe-alignment=/gpfs/scratch/bsc59/bsc59829/Spatial/COVID/Metadata/HRC/Images/Control2.json
