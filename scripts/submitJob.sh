#!/bin/bash
#
# submitJob.sh: Job submission script for CCA-based scRNA-Seq pipeline
# By Kristin Muench
# 2019.07.20
#
# SLURM submission flags:
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=124G
#SBATCH --partition=interactive
#SBATCH --job-name=idealPipelineAllRes
#SBATCH --time=6-00:00:00
#SBATCH --mail-user=kmuench@stanford.edu
#SBATCH --mail-type=ALL
#
# Many filepaths (e.g. $METADAT_16p_SC_MS) are located in my .bash_profile
#
# Paper motivating Seurat3: https://www-cell-com.stanford.idm.oclc.org/cell/fulltext/S0092-8674(19)30559-8
#
# Who did extraction metadata: /labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190720_whoDidExtraction.csv
#
# Pipeline done according to how I think it should be done: $OUTPUT_16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500
# Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/20190716_remakePipeline
# No filtering of data: Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_normalization_noFilter
# No normalizing of data: Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_noNormalizing_filter1000filter7500
# just smple H: Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_H
# jsut sample I: Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_I
#
# Example run:
# qsub submitJob.sh
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# prepare modules
module purge
module load r/3.5.0
module load miniconda/3 # to get access to UMAP

# load in paths
#barcodeSampleLUT=/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190326_troubleshootRunthrough/determineCellSex/allNewID_cells.csv # could be demultiplexed IDs or not
outputDir="$OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_noNormalizing_filter1000filter7500" # directory where output stored

# declare for documentation purposes what these variables are
echo 'Barcode Sample Lookup Table Location: ' $barcodeSampleLUT
echo 'Output Directory: ' $outputDir

# Import data
## Order of variables: (1)output loc , (2) Metadata.csv , (3) Folder containin .mtx data from 10x  # (4) (if applicable) barcodeSampleLUT     e.g. $barcodeSampleLUT
echo 'Making variables...'
Rscript 1_initializeVars.R $outputDir $METADAT_16p_SC_MS $MTX_16p_SC_MS  

# Find clusters and visualize batch effects
## (3) path To Data Combined
echo 'Finding clusters...'
Rscript 2_cluster.R $outputDir $METADAT_16p_SC_MS "$outputDir/initializeVars/seuratObj_data.combined.RData"







# Merge the individual samples and perform CCA, visualize output to decide on number of CCs to use
## Order of variables: (1)Folder containing input from makeVars.R, (2) output loc   
#outputDir_makeVars="$outputDir/makeVars/allVars.RData"
#echo 'Calculating CCA...'
#echo 'Importing this makeVars.RData file:' $outputDir_makeVars
#Rscript calcCCA.R $outputDir_makeVars $outputDir

# Align based on CCs, do tSNE, find clusters
## Order of variables: 
### (1) CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#echo 'Calculating Clusters...'
#Rscript clusters.R "$outputDir/calcCCA/data.combined_multiCCA_metadata.RData" $outputDir 40 1.2

# Visualize output
## Order of variables: 
### (1) File containing CCA-bearing variable set up in cca_calcCCA.R, 
### (2) output container folder,
### (3) chosen number of CCs
### (4) the chosen resolution
#echo 'Making useful visualizations...'
#outputDir_clusters="$outputDir/clusters/data.combined_withClust_r1.2_CC40.RData"
#echo 'Importing this RData file with data and clusters:'
#Rscript groupCompare.R $outputDir_clusters $outputDir 40 1.2
