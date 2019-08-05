# 3_analyzeClustersUsingLabels.R
# Kristin Muench
# 2019.08.-2
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


# !!! FIX UP TO IMPORT DATA FROM STEP 3 AND THEN VISUALIZE VLUSTER LABELS - GOOD PLACE FOR PIE CHART?

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
allVarsPath <- args[3]
singleCellClusterMarkersPath <- args[4]
demuxId <- args[5]

# Import data
#allVarsPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/from20190720_excitCluster/cluster/allVars_res1.2.RData'
load(allVarsPath)

#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
metadata<-read.csv(metadataPath)

#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/from20190720_excitCluster'

# singleCellClusterMarkersPath <- '/labs/tpalmer/resources/singleCellClusterMarkers/20190729_markerGenes_comprehensive.csv'
singleCellClusterMarkers <- read.csv(singleCellClusterMarkersPath)


# Create folder to store output
setwd(outputDir)
subDir <- 'analyzeClustersUsingLabels'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}


# Import demultiplexed IDs if possible
argsLen <- length(args);

outfile <- if (argsLen < 5) {
  print('Using non-demultiplexed IDs...')
}else{
  demuxID <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/sexDemux/mapCellsToNewIDs.RData'
  print(paste0('Using demultiplexed IDs from from ', demuxID))
  load(demuxID)
  
  data.combined@meta.data$SampleID <- mapCellsToNewIDs[match(row.names(data.combined@meta.data), mapCellsToNewIDs$cells),'demuxID']
  
  # Generate some UMAPs on the basis of sex
  bc_sexDemuxFolder <- 'batchEffects_sexDemux'
  dir.create(file.path(outputDir, subDir, bc_sexDemuxFolder))
  setwd(file.path(outputDir, subDir, bc_sexDemuxFolder))
  
  p1 <- DimPlot(data.combined, reduction = "umap", group.by = 'SampleID')
  p_split <- DimPlot(data.combined, reduction = "umap", split.by = 'SampleID', ncol=3)
    
  ggsave(filename = paste0('umap_batch_demuxID.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 20, height=20, units ='cm')
  ggsave(filename = paste0('umap_batch_splitby_demuxID.pdf'), plot = p_split, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 40, height=40, units ='cm')
  
}


# Dot plot to visualize major markers 
data.combined@active.ident <- factor(data.combined@meta.data$seurat_clusters, 
                                     levels = c('1','0','25','16','19',
                                                '3','21','7','4','22',
                                                '23','6','9','14','15', 
                                                '20','13','17','11','8',
                                                '5','2','10','12','18') )# create ordering that will be pretty
markers.to.plot <- c('Tuba1b', 'Pax6', 'Sox2', 'Top2a', 'Slbp', 'Rrm2', 'Eomes', 'Neurog2', 'Neurod2', 'Mapt', 'Bcl11b', 'Crym', 'Ldb2','Reln')
dp <- DotPlot(data.combined, features = rev(markers.to.plot), dot.scale = 8) + RotatedAxis()

ggsave(filename = paste0('visMarkersForClusters.pdf'), 
       plot = dp, device='pdf', 
       path = file.path(outputDir, subDir), width = 20, height=30, units ='cm')

# SCRATCH ## Version of this with group.by CATEGORY!!!!!
data.combined@meta.data$GeneralCategory <- 'Radial Glia'
data.combined@meta.data[which(data.combined@meta.data$seurat_clusters %in% c()),'GeneralCategory'] <- 'Intermediate Progenitors'







# # BONUS - Generate list of cells to use in subset analysis from imported cluster annotation

pathToClusterLabels <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/cluster/clusterMarkers_res1.2/20190729_clusterAssignments_fullDataset_res1.2.csv'
clusterLabels <- read.csv(pathToClusterLabels)
clustersToUse <- clusterLabels[which(clusterLabels$ExcitatoryLineage == TRUE),'Cluster']

cellsToUse <- row.names(data.combined@meta.data[which(data.combined@meta.data$seurat_clusters %in% clustersToUse),])
setwd(file.path(outputDir, subDir))
save(cellsToUse, file = 'cells_excitatorySubsetAnalysis.RData')



#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save(ls(pattern='markersPresent_*'), file = 'markersPresent_all.RData')
save.image(file = paste0(paste0("allVars.RData")))
save(data.combined, file = 'seuratObj_data.combined.RData')

print('~*~ All done! ~*~')









