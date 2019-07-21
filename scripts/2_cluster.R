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

print(paste0('Seurat package Version: ', packageVersion('Seurat')) )
print(paste0('Cowplot package Version: ', packageVersion('cowplot')) )
print(paste0('Ggplot2 package Version: ', packageVersion('ggplot2')) )

# import from command line
args <- commandArgs(TRUE)
outputDir <- args[1]
metadataPath <- args[2]
dataCombinedPath <- args[3]
resToTry <- args[4]

# Import data
#dataCombinedPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/s1s2_I/initializeVars/seuratObj_data.combined.RData'
load(dataCombinedPath)

#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
metadata<-read.csv(metadataPath)

#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/s1s2_I/'


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
myMetadata <- merge(metadata[,metadataCols], myCells, by.x='SampleID', by.y='Identity')
rownames(myMetadata) <- myMetadata$cellNames
data.combined <- AddMetaData(data.combined, myMetadata, col.name = metadataCols) ## add metadata to Seurat object


# Perform clustering

## Declare that this is integrated analysis
DefaultAssay(data.combined) <- "integrated" # indicate that this is an integrated dataset

## Run the standard workflow for visualization and clustering
print('Scaling data, running PCA...')
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 70, verbose = FALSE)

## t-SNE and Clustering
print('Running UMAP and finding neighbors...')
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:49)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:49)

save(data.combined, file = 'seuratObj_data.combined_beforeClustering.RData')

## Find clusters for list of resolutions
print('Finding clusters...')

for (r in resToTry){
  print(paste0('Finding clusters for resolution = ', r))
  data.combined <- FindClusters(data.combined, resolution = r)
  
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
  
  
  # Identify cell type markers
  ## Generate lists and save
  numClust <- max(as.numeric(data.combined@meta.data[,paste0('integrated_snn_res.',r)]))
  print(paste0('The numbers of clusters found is ', numClust))
  for (c in c(0:(numClust-1))){
    print(paste0('Finding markers for cluster ', c, '...'))
    
    tmp_markers <- FindConservedMarkers(data.combined, ident.1 = as.numeric(c), grouping.var = "cond", verbose = FALSE)
    
    write.csv(tmp_markers, file = paste0('clustMarker_',c, '.csv'))
    assign(paste0('clustMarker_', c) , tmp_markers)
    
    visGenes <- FeaturePlot(data.combined, features = row.names(head(tmp_markers)), min.cutoff = "q9", path = file.path(outputDir, subDir, subDir_clusterMarkers), width = 25, height=20, units ='cm')
    ggsave(paste0('featurePlot_', c, '.pdf'), plot = visGenes, device = 'pdf')
  }
  
  #objectsToSave <- ls()[grep('clustMarker_',ls())]
  #save( eval(parse(text=objectsToSave)) , file = paste0('clusterMarkers_res', r,'.RData'))
  
  # Visualize batch effects
  print(paste0('Visualizing batch effects for resolution ', r, '...'))
  
  ## Create subdirectory to store data
  subDir_batchEffects <- paste0('batchEffects_res', r)
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
    
    ggsave(filename = paste0('umap_batch_', b,'.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 20, height=20, units ='cm')
    ggsave(filename = paste0('umap_batch_splitby_', b,'.pdf'), plot = p_split, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 40, height=20, units ='cm')
  }
  
  
  
}



#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save.image(file = paste0("allVars.RData"))
save(data.combined, file = 'seuratObj_data.combined.RData')


# Dot plot to visualize major markers
markers.to.plot <- c('Tuba1b', 'Pax6', 'Sox2', 'Top2a', 'Slbp', 'Rrm2', 'Eomes', 'Neurog2', 'Neurod2', 'Mapt', 'Bcl11b', 'Crym', 'Ldb2','Reln', 'Gad2', 'Sst', 'Trem2', 'Igfbp7', 'Otx2', 'Tcf7l2', 'Col1a2', 'Lum', 'Hbb.bh1')
dp <- DotPlot(data.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "cond") + RotatedAxis()

ggsave(filename = paste0('visMarkersForClusters_res', resToTry,'.pdf'), plot = dp, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 20, height=40, units ='cm')

print('~*~ All done! ~*~')