# 2_cluster.R
# Kristin Muench
# 2019.07.20
# Pipeline to run Seurat to analyze Mouse scRNA-Seq data from a 10x Chromium.
# This is the first step in the pipeline.
# Seurat v3.0.2 is used
#
# 1. InitializeVars - Import data from 10X, Make combined object
# 2. Cluster - Scale, run PCA, run UMAP, Find Neighbors, Find Clusters, Visualize Batch Effects, Identify Cell Type Markers
# 3. ComparePopCounts - compare population sizes between factors of interest
# 4. CompareDE - compare differential expression between factors of interest
#
# S. Script to demultiplex on the basis of sex, and identify putative doublets using co-expression of Xist and Ymarker genes
# C. Script for comparing cluster stability when clusters made using different methods?
#
# Explains the basic QC steps: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
#
# Following tutorial here: https://satijalab.org/seurat/v3.0/immune_alignment.html
# 
# Integration explained: https://satijalab.org/seurat/v3.0/integration.html
#
# # # # # # # # # # # # # # # # # #

# Load needed library
library('Seurat') 
library('cowplot')
library('ggplot2')

packageVersion('Seurat')
packageVersion('cowplot')
packageVersion('ggplot2')

# import from command line
args <- commandArgs(TRUE)
outputDir <- args[1]
metadataPath <- args[2]
dataCombinedPath <- args[3]

# Import data
#dataCombinedPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190716_remakePipeline/initializeVars/seuratObj_data.combined.RData'
load(dataCombinedPath)

#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
metadata<-read.csv(metadataPath)

#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190716_remakePipeline/'


# Create folder to store output
setwd(outputDir)
subDir <- 'cluster'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}


# add metadata to Seurat object
myCells <- data.frame(cellNames = row.names(data.combined@meta.data), Identity = data.combined@meta.data$orig.ident)
metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter', 'OrderOfLitterExtraction')
myMetadata <- merge(myCells, metadata[,metadataCols], by.x='Identity', by.y='SampleID')
rownames(myMetadata) <- myMetadata$cellNames
data.combined <- AddMetaData(data.combined, myMetadata[,-c(1:2)], col.name = metadataCols[-c(1:2)]) ## add metadata to Seurat object


# Perform clustering

## Declare that this is integrated analysis
DefaultAssay(data.combined) <- "integrated" # indicate that this is an integrated dataset

## Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 70, verbose = FALSE)

## t-SNE and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:40)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:40)


## Find clusters for list of resolutions

resToTry <- c(0.8, 1.2, 2, 7, 10)

for (r in resToTry){
  print(paste0('Finding clusters for resolution = ', resToTry))
  data.combined <- FindClusters(data.combined, resolution = 1.2)
  
  ##Visualize
  p_labels <- DimPlot(data.combined, reduction = "umap", label = TRUE)
  ggsave(filename = paste0(paste0('umap_all_labels_res', r, '.pdf')), plot = p_labels, device='pdf', path = file.path(outputDir, subDir), width = 20, height=20, units ='cm')
  
  ## Create subdirectory to store data
  subDir_clusterMarkers <- paste0('clusterMarkers_res', r)
  setwd(file.path(outputDir, subDir))
  if (file.exists(subDir_clusterMarkers)){
    setwd(file.path(outputDir, subDir,subDir_clusterMarkers))
  } else {
    dir.create(file.path(outputDir, subDir, subDir_clusterMarkers))
    setwd(file.path(outputDir, subDir, subDir_clusterMarkers))
  }
  
  ## Generate lists and save
  for (c in max(data.combined@meta.data$res1.2)){
    print(paste0('Finding markers for cluster ', c, '...'))
    
    tmp_markers <- FindConservedMarkers(immune.combined, ident.1 = c, grouping.var = "stim", verbose = FALSE)
    head(tmp_markers)
    
    write.csv(tmp_markers)
    assign(paste0('clustMarker_', c) , tmp_markers)
  }
  
  
}


# Visualize batch effects

## Create subdirectory to store data
subDir_batchEffects <- 'batchEffects'
setwd(file.path(outputDir, subDir))
if (file.exists(subDir_batchEffects)){
  setwd(file.path(outputDir, subDir,subDir_batchEffects))
} else {
  dir.create(file.path(outputDir, subDir, subDir_batchEffects))
  setwd(file.path(outputDir, subDir, subDir_batchEffects))
}

## Generate visualizations
for (b in metadataCols){
  p1 <- DimPlot(data.combined, reduction = "umap", group.by = b)
  p_split <- DimPlot(data.combined, reduction = "umap", split.by = b)
  
  ggsave(filename = paste0('umap_batch_', b,'.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir), width = 20, height=20, units ='cm')
  ggsave(filename = paste0('umap_batch_splitby_', b,'.pdf'), plot = p_split, device='pdf', path = file.path(outputDir, subDir), width = 20, height=20, units ='cm')
}


# Identify cell type markers





##  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save.image(file = paste0("allVars.RData"))
save(data.combined, file = 'seuratObj_data.combined.RData')
save(ls()[grep('clustMarker_',ls())], file = 'clusterMarkers.RData')
