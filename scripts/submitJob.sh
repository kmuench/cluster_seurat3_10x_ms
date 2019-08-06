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
#SBATCH --job-name=e_sc
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
# RC normalization: Test pipeline: $OUTPUT_16p/201907_cluster_seurat_10x_ms/s1s2_RCnorm
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
outputDir_all="$OUTPUT_16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500" # directory where output stored

# declare for documentation purposes what these variables are
echo 'Barcode Sample Lookup Table Location: ' $barcodeSampleLUT
echo 'Output Directory: ' $outputDir_all

# Import data
## Order of variables: (1)output loc , (2) Metadata.csv , (3) Folder containin .mtx data from 10x  # (4) (if applicable) barcodeSampleLUT     e.g. $barcodeSampleLUT
#echo 'Making variables...'
#Rscript 1_initializeVars.R $outputDir_all $METADAT_16p_SC_MS $MTX_16p_SC_MS  

# Find clusters and visualize batch effects
## (3) path To Data Combined, (4) resolution param
## To change in the future: pass in resolution so you can run different resolutions as different jobs (in parallel)
#echo 'Finding clusters...'
#Rscript 2_cluster.R $outputDir_all $METADAT_16p_SC_MS "$outputDir_all/initializeVars/seuratObj_data.combined.RData" 1.2

# Find clusters and visualize batch effects
## (3) path To Data Combined, (4) resolution param
## To change in the future: pass in resolution so you can run different resolutions as different jobs (in parallel)
#echo 'Demultiplexing sex...'
#Rscript S_sexDemux.R $outputDir_all $METADAT_16p_SC_MS "$outputDir_all/cluster/seuratObj_data.combined.RData"

# DE Analysis
#Rscript 5_differentialExpression.R $outputDir_all $METADAT_16p_SC_MS "$outputDir_all/analyzeClusters/seuratObj_data.combined.RData" "$outputDir_all/cluster/clusterMarkers_res1.2/20190729_clusterAssignments_fullDataset_res1.2.csv"  



# BEGIN THE PART OF THE PIPELINE THAT IS SUBSETTING EXCITATORY CELLS #
outputDir_excit="$OUTPUT_16p/201907_cluster_seurat_10x_ms/from20190720_excitCluster" # directory where output stored

# Find clusters and visualize batch effects
## (3) path To Data Combined, (4) resolution param
echo 'Finding clusters...'
#Rscript 2_cluster.R $outputDir_excit $METADAT_16p_SC_MS "$outputDir_all/sexDemux/seuratObj_data.combined_demux.RData" 1.2 "$outputDir_all/analyzeClusters/cells_excitatorySubsetAnalysis.RData"

# Analyze the clusters generally
## 
#Rscript 3_analyzeClusters.R $outputDir_excit $METADAT_16p_SC_MS "$outputDir_excit/cluster/seuratObj_data.combined_res1.2.RData" "/labs/tpalmer/resources/singleCellClusterMarkers/20190729_markerGenes_comprehensive.csv" "$outputDir_all/sexDemux/mapCellsToNewIDs.RData"

# DE analysis
Rscript 5_differentialExpression.R $outputDir_excit $METADAT_16p_SC_MS "$outputDir_excit/analyzeClusters/seuratObj_data.combined.RData" "$outputDir_excit/cluster/20190802_clusterAssignments_excitOnly_res1.2.csv" 

echo 'All done with job script!'
