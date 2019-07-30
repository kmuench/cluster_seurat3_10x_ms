# 3_analyzeClusters.R
# Kristin Muench
# 2019.07.29
# Pipeline to run Seurat to analyze Mouse scRNA-Seq data from a 10x Chromium.
# This is the step in which you visualize what cell types are.
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
singleCellClusterMarkersPath <- args[4]

# Import data
#dataCombinedPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/cluster/allVars_res1.2.RData' # allVars before I rerun to make a res-specific allVars
load(dataCombinedPath)

#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
metadata<-read.csv(metadataPath)

#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500'

# singleCellClusterMarkersPath <- '/labs/tpalmer/resources/singleCellClusterMarkers/20190729_markerGenes_comprehensive.csv'
singleCellClusterMarkers <- read.csv(singleCellClusterMarkersPath)

# Create folder to store output
setwd(outputDir)
subDir <- 'analyzeClusters'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}


# Dot plot to visualize major markers
markers.to.plot <- c('Tuba1b', 'Pax6', 'Sox2', 'Top2a', 'Slbp', 'Rrm2', 'Eomes', 'Neurog2', 'Neurod2', 'Mapt', 'Bcl11b', 'Crym', 'Ldb2','Reln', 'Gad2', 'Sst', 'Trem2', 'Igfbp7', 'Otx2', 'Tcf7l2', 'Col1a2', 'Lum', 'Hbb-bh1')
dp <- DotPlot(data.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
              split.by = "cond") + RotatedAxis()

ggsave(filename = paste0('visMarkersForClusters_res', resToTry,'.pdf'), plot = dp, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 20, height=40, units ='cm')


# For each marker list, are there any known cluster markers in the top 50 positive FC genes?

## create subfolder to hold all files
subDir_countMarkers <- 'countMarkers'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir, subDir_countMarkers))
} else {
  dir.create(file.path(outputDir, subDir, subDir_countMarkers))
  setwd(file.path(outputDir, subDir, subDir_countMarkers))
}

## function to identify whether or not known markers in the list
countMarkers <- function(clust, singleCellClusterMarkers){
  # what genes are the top in this cluster?
  topMarkers <- row.names(clust[rev(order(clust$WT.SAL_avg_logFC))[1:50],])
  whichMarkersPresent <- singleCellClusterMarkers[which(singleCellClusterMarkers$Transcript %in% topMarkers),]
  
  numMarkersPerCategory <- apply(singleCellClusterMarkers[,!(names(singleCellClusterMarkers) %in% c('Transcript','Citations', 'Other.Notes', 'Microglia') )],
                                 2, FUN = function(x) sum(x, na.rm=TRUE))
  
  # score - how many per category?
  scores <- apply(whichMarkersPresent[,!(names(whichMarkersPresent) %in% c('Transcript','Citations', 'Other.Notes', 'Microglia') )],
                  2, FUN = function(x) sum(x, na.rm=TRUE))
  
  scores_full <- merge(scores, numMarkersPerCategory, by='row.names')
  colnames(scores_full) <- c('category', 'count', 'totalMarkersInCategory')
  scores_full$percents <- (scores_full$count/scores_full$totalMarkersInCategory)*100
  
  scores_full <- scores_full[rev(order(scores_full$percents)),]
  
  print(scores_full[c(1:10),])

  # save csv for documentation purposes
  write.csv(file = paste0('whichMarkersPresent_', deparse(substitute(clust)), '.csv'), whichMarkersPresent)
  write.csv(file = paste0('numberOfHitsPerCategory_', deparse(substitute(clust)), '.csv'), scores)
  
  # return the list if desired
  return( whichMarkersPresent )
}

## deploy on all marker genes
markerGeneLists <- ls(pattern='clustMarker*')

for (m in markerGeneLists) {
  print(paste0('Marker scores for cluster ', m) )
  assign(paste0('markersPresent_', m), countMarkers(eval(parse(text = m)), singleCellClusterMarkers) )
}




#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save(ls(pattern='markersPresent_*'), file = 'markersPresent_all.RData')
save.image(file = paste0(paste0("allVars.RData")))
save(data.combined, file = 'seuratObj_data.combined.RData')

print('~*~ All done! ~*~')

# # BONUS - Generate list of cells to use in subset analysis from imported cluster annotation

pathToClusterLabels <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/cluster/clusterMarkers_res1.2/20190729_clusterAssignments_fullDataset_res1.2.csv'
clusterLabels <- read.csv(pathToClusterLabels)
clustersToUse <- clusterLabels[which(clusterLabels$ExcitatoryLineage == TRUE),'Cluster']

cellsToUse <- row.names(data.combined@meta.data[which(data.combined@meta.data$seurat_clusters %in% clustersToUse),])
setwd(file.path(outputDir, subDir))
save(cellsToUse, file = 'cells_excitatorySubsetAnalysis.RData')










#  ## SCRATCH
# # add metadata to Seurat object
# myCells <- data.frame(cellNames = row.names(data.combined@meta.data), Identity = data.combined@meta.data$orig.ident)
# metadataCols <- c('SampleID','CellsPerSample','SurgeryDate','Condition', 'Genotype', 'Litter', 'OrderOfLitterExtraction')
# myMetadata <- merge(metadata[,metadataCols], myCells, by.x='SampleID', by.y='Identity')
# rownames(myMetadata) <- myMetadata$cellNames
# data.combined <- AddMetaData(data.combined, myMetadata, col.name = metadataCols) ## add metadata to Seurat object
# 
# 
# # Perform clustering
# 
# ## Declare that this is integrated analysis
# DefaultAssay(data.combined) <- "integrated" # indicate that this is an integrated dataset
# 
# ## Run the standard workflow for visualization and clustering
# print('Scaling data, running PCA...')
# data.combined <- ScaleData(data.combined, verbose = FALSE)
# data.combined <- RunPCA(data.combined, npcs = 70, verbose = FALSE)
# 
# ## t-SNE and Clustering
# print('Running UMAP and finding neighbors...')
# data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:49)
# data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:49)
# 
# save(data.combined, file = paste0('seuratObj_data.combined_beforeClustering_res', resToTry, '.RData') )
# 
# ## Find clusters for list of resolutions
# print('Finding clusters...')
# 
# for (r in resToTry){
#   print(paste0('Finding clusters for resolution = ', r))
#   data.combined <- FindClusters(data.combined, resolution = as.numeric(r))
#   
#   ##Visualize
#   p_labels <- DimPlot(data.combined, reduction = "umap", label = TRUE)
#   ggsave(filename = paste0(paste0('umap_all_labels_res', r, '.pdf')), plot = p_labels, device='pdf', path = file.path(outputDir, subDir), width = 20, height=20, units ='cm')
#   
#   ## Create subdirectory to store data
#   subDir_clusterMarkers <- paste0('clusterMarkers_res', r)
#   setwd(file.path(outputDir, subDir))
#   if (file.exists(subDir_clusterMarkers)){
#     setwd(file.path(outputDir, subDir,subDir_clusterMarkers))
#   } else {
#     dir.create(file.path(outputDir, subDir, subDir_clusterMarkers))
#     setwd(file.path(outputDir, subDir, subDir_clusterMarkers))
#   }
#   
#   
#   # Identify cell type markers
#   ## Generate lists and save
#   numClust <- max(as.numeric(data.combined@meta.data[,paste0('integrated_snn_res.',r)]))
#   print(paste0('The numbers of clusters found is ', numClust))
#   for (c in c(0:(numClust-1))){
#     print(paste0('Finding markers for cluster ', c, '...'))
#     
#     tmp_markers <- FindConservedMarkers(data.combined, ident.1 = as.numeric(c), grouping.var = "cond", verbose = FALSE)
#     
#     write.csv(tmp_markers, file = paste0('clustMarker_',c, '.csv'))
#     assign(paste0('clustMarker_', c) , tmp_markers)
#     
#     visGenes <- FeaturePlot(data.combined , features = row.names(head(tmp_markers)) , min.cutoff = "q9")
#     ggsave(paste0('featurePlot_', c, '.pdf'), plot = visGenes, device = 'pdf', path = file.path(outputDir, subDir, subDir_clusterMarkers), width = 25, height=20, units ='cm')
#   }
#   
#   #objectsToSave <- ls()[grep('clustMarker_',ls())]
#   #save( eval(parse(text=objectsToSave)) , file = paste0('clusterMarkers_res', r,'.RData'))
#   
#   # Visualize batch effects
#   print(paste0('Visualizing batch effects for resolution ', r, '...'))
#   
#   ## Create subdirectory to store data
#   subDir_batchEffects <- paste0('batchEffects_res', r)
#   setwd(file.path(outputDir, subDir))
#   if (file.exists(subDir_batchEffects)){
#     setwd(file.path(outputDir, subDir,subDir_batchEffects))
#   } else {
#     dir.create(file.path(outputDir, subDir, subDir_batchEffects))
#     setwd(file.path(outputDir, subDir, subDir_batchEffects))
#   }
#   
#   ## Generate visualizations
#   for (b in metadataCols){
#     p1 <- DimPlot(data.combined, reduction = "umap", group.by = b)
#     p_split <- DimPlot(data.combined, reduction = "umap", split.by = b)
#     
#     ggsave(filename = paste0('umap_batch_', b,'.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 20, height=20, units ='cm')
#     ggsave(filename = paste0('umap_batch_splitby_', b,'.pdf'), plot = p_split, device='pdf', path = file.path(outputDir, subDir, subDir_batchEffects), width = 40, height=20, units ='cm')
#   }
#   
#   
#   
# }
# 
# 
# 
# 
# 
# 
# 
