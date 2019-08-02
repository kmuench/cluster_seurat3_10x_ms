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
library('dplyr')

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
subDir <- 'analyzeClusters'

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
  p1_split <- DimPlot(data.combined, reduction = "umap", split.by = 'SampleID', ncol=3)
  
  p2 <- DimPlot(data.combined, reduction = "umap", group.by = 'Sex')
  p2_split <- DimPlot(data.combined, reduction = "umap", split.by = 'Sex', ncol=3)
    
  ggsave(filename = paste0('umap_batch_demuxID.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 20, height=20, units ='cm')
  ggsave(filename = paste0('umap_batch_splitby_demuxID.pdf'), plot = p1_split, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 40, height=40, units ='cm')
  
  ggsave(filename = paste0('umap_batch_Sex.pdf'), plot = p2, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 20, height=20, units ='cm')
  ggsave(filename = paste0('umap_batch_splitby_Sex.pdf'), plot = p2_split, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 40, height=40, units ='cm')
  
}

# # SPECIAL NOTE AUG 2 2019
# # I decided to take out Cluster 26 from the "from20190720_excitCluster" analysis-uncomment to remove that cluster 
# 
# data.combined_full <- data.combined
# idents <- levels(data.combined_full@active.ident)
# data.combined <- subset(data.combined, idents = idents[!(idents %in% c('26','24'))] ) # also see ident.remove



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
countMarkers <- function(clust, singleCellClusterMarkers, m){
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
  write.csv(file = paste0('whichMarkersPresent_', m, '.csv'), whichMarkersPresent)
  write.csv(file = paste0('numberOfHitsPerCategory_', m, '.csv'), scores)
  
  # return the list if desired
  return( whichMarkersPresent )
}

## deploy on all marker genes
markerGeneLists <- ls(pattern='clustMarker*')

for (m in markerGeneLists) {
  print(paste0('Marker scores for cluster ', m) )
  assign(paste0('markersPresent_', m), countMarkers(eval(parse(text = m)), singleCellClusterMarkers, m) )
}



# Are there any population size differences between groups?

## create subfolder to hold all files
subDir_popTallies <- 'popTallies'

if (file.exists(subDir_popTallies)){
  setwd(file.path(outputDir, subDir, subDir_popTallies))
} else {
  dir.create(file.path(outputDir, subDir, subDir_popTallies))
  setwd(file.path(outputDir, subDir, subDir_popTallies))
}

## perform tallying
tallyClust <- data.frame(data.combined@meta.data %>% group_by(SampleID, seurat_clusters) %>% tally())
totalCellsPerSample <- data.frame(data.combined@meta.data %>% group_by(SampleID) %>% count())
tallyClust$total <- data.frame(totalCellsPerSample[match(tallyClust$SampleID, totalCellsPerSample$SampleID),'n'])[,1]
tallyClust$percent <- (tallyClust$n / tallyClust$total)*100
neutralColors1 <- c('#009b76', '#eaab00') 
neutralColors2 <- c('#009b76', '#eaab00') 
tableauColors <- c('#4e79a7', '#59a14f', '#9c755f', '#f28e2b', '#edc948', '#bab0ac', '#e15759', '#b07aa1', '#76b7b2', '#ff9da7')

## ensure levels are what you want
data.combined@meta.data$cond <- factor(data.combined@meta.data$cond, levels = c('WT.SAL', 'WT.LPS', 'HET.SAL', 'HET.LPS'))
data.combined@meta.data$Sex <- factor(data.combined@meta.data$Sex, levels = c('M', 'F'))
data.combined@meta.data$Condition <- factor(data.combined@meta.data$Condition, levels = c('SAL', 'LPS'))
data.combined@meta.data$Genotype <- factor(data.combined@meta.data$Genotype, levels = c('WT', 'HET'))

                                                                              

barPlotComparePop_1fac <- function(tallyClust, metaData, chosenFactor, chosenClust, neutralColors){
  # add factor to list
  tallyClust$factor <- metaData[match(tallyClust$SampleID, metaData$SampleID), chosenFactor]
  
  # just pull out data for the cluster of interest
  tallyClust_subset <- tallyClust %>% filter(seurat_clusters %in% chosenClust)
  
  # one-way anova
  res.aov <- aov(percent ~ factor, data = tallyClust_subset)
  print(summary(res.aov))
  write.table(unclass(summary(res.aov)), file = paste0('stats_group',chosenFactor, '_aov', paste(chosenClust, sep="", collapse="_"), '.csv'))
  
  # generate mean and sd
  plotData <- tallyClust_subset %>% group_by(factor) %>% summarize(myMean = mean(percent), mySd = sd(percent))
  
  p <- ggplot(plotData, aes(x=factor, y=myMean, fill=factor) ) +
    geom_bar(position=position_dodge(), stat="identity",
             color = "black", size = 0.3) +
    geom_errorbar(aes(ymin=myMean-mySd, ymax=myMean+mySd), width=.2,
                  position=position_dodge(.9)) +
    theme_bw() +
    xlab(chosenFactor) +
    ylab("Percent of Cells in Category")+
    ggtitle(paste0('Percent of Cells in Excitatory Lineage in Cluster ', paste(chosenClust, sep="", collapse="_"))) +
    scale_fill_manual(values=neutralColors )
  
  ggsave(filename = paste0('barPlot_1factor', chosenFactor, '_clust', paste(chosenClust, sep="", collapse="_"), '.pdf'), plot = p, device='pdf',
         width = 20, height=20, units ='cm')
  
  return(res.aov)
}

barPlotComparePop_2fac <- function(tallyClust, metaData, chosenFactor_A, chosenFactor_B, chosenClust, neutralColors){
  # add factor to list
  tallyClust$factor_A <- metaData[match(tallyClust$SampleID, metaData$SampleID), chosenFactor_A]
  tallyClust$factor_B <- metaData[match(tallyClust$SampleID, metaData$SampleID), chosenFactor_B]
  
  # just pull out data for the cluster of interest
  tallyClust_subset <- tallyClust %>% filter(seurat_clusters %in% chosenClust)
  
  # one-way anova
  res.aov <- aov(percent ~ factor_A * factor_B, data = tallyClust_subset)
  print(summary(res.aov))
  write.table(unclass(summary(res.aov)), file = paste0('stats_',chosenFactor_A, '_', chosenFactor_B, '_aov', paste(chosenClust, sep="", collapse="_"), '.csv'))
  
  # generate mean and sd
  plotData <- tallyClust_subset %>% group_by(factor_A, factor_B) %>% summarize(myMean = mean(percent), mySd = sd(percent))
  
  p <- ggplot(plotData, aes(x=factor_B, y=myMean, fill=factor_A) ) +
    geom_bar(position=position_dodge(), stat="identity",
             color = "black", size = 0.3) +
    geom_errorbar(aes(ymin=myMean-mySd, ymax=myMean+mySd), width=.2,
                  position=position_dodge(.9)) +
    theme_bw() +
    xlab(chosenFactor_A) +
    ylab("Percent of Cells in Category")+
    ggtitle(paste0('Percent of Cells in Excitatory Lineage in Cluster ', paste(chosenClust, sep="", collapse="_"))) +
    scale_fill_manual(values= neutralColors )
  
  ggsave(filename = paste0('barPlot_2factor', chosenFactor_A, '_', chosenFactor_B, '_clust', paste(chosenClust, sep="", collapse="_"), '.pdf'), plot = p, device='pdf',
         width = 30, height=20, units ='cm')
  
  return(res.aov)
}

for (clust in c(levels(tallyClust$seurat_clusters)) ){
  print(paste0('Analysis for ', clust))
  # create folder
  dirName <- paste0('cluster', clust)
  if (file.exists(file.path(outputDir, subDir, subDir_popTallies, dirName))){
    setwd(file.path(outputDir, subDir, subDir_popTallies, dirName))
  } else {
    dir.create(file.path(outputDir, subDir, subDir_popTallies, dirName))
    setwd(file.path(outputDir, subDir, subDir_popTallies, dirName))
  }
  print(getwd())
  
  # one factor: factor is...
  ## cond
  print('One-factor analysis for condition...')
  barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'cond', clust, tableauColors[c(1:4)])

  ## sex
  print('One-factor analysis for sex...')
  barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', clust, tableauColors[c(6,5)])
  
  ## Condition
  print('One-factor analysis for Condition...')
  barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Condition', clust, tableauColors[c(1,2)])
  
  ## Genotype
  print('One-factor analysis for Genotype...')
  barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Genotype', clust, tableauColors[c(7:8)])
  
  
  # two factor: factors are...
  ## Inject * Geno
  barPlotComparePop_2fac(tallyClust, metaData, 'Genotype', 'Condition', chosenClust, tableauColors[c(7:8)])
  
  ## Cond * Sex
  barPlotComparePop_2fac(tallyClust, metaData, 'Sex', 'cond', chosenClust, tableauColors[c(6,5)])
  
}

# Follow up on groups of cells

## Do F have different transits through IPC stage?

# MvF in  11/23/22? 
setwd(file.path(outputDir, subDir, subDir_popTallies))
data.combined@meta.data$IPC_Path_A <- FALSE
data.combined@meta.data[which(data.combined@meta.data$seurat_clusters %in% c('11', '23', '22')),'IPC_Path_A'] <- TRUE

barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('11', '23', '22'), tableauColors[c(6,5)])
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('23', '22'), tableauColors[c(6,5)])


# MvF in 9/14/15/20?
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('9', '14', '15', '20'), tableauColors[c(6,5)])
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('14', '15', '20'), tableauColors[c(6,5)])

#MvF in 9/13
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('9', '13'), tableauColors[c(6,5)])

# This just looks like a group
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'Sex', c('11', '23', '13'), tableauColors[c(6,5)])


# Does condition affect postmitotic neurons
  # 12, 18
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'cond', c('12', '18'), tableauColors)

  # 10, 12, 18
barPlotComparePop_1fac(tallyClust, data.combined@meta.data, 'cond', c('10', '12', '18'), tableauColors)


#What's the difference between 9 and 13?
diff_between_9_and_13 <- FindMarkers(data.combined, ident.1 = "9", ident.2 = "13")
write.csv(diff_between_9_and_13, file ='diff_between_9_and_13.csv')

# Just visualize clusters
## Visualize 9:
# DimPlot(data.combined, cells = row.names(data.combined@meta.data[which(data.combined@meta.data$seurat_clusters == '9'),] ) )


#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save(ls(pattern='markersPresent_*'), file = 'markersPresent_all.RData')
save.image(file = paste0(paste0("allVars.RData")))
save(data.combined, file = 'seuratObj_data.combined.RData')

print('~*~ All done! ~*~')






