# 5_differentialExpression.R
# Kristin Muench
# 2019.08.05
# Pipeline to run Seurat to analyze Mouse scRNA-Seq data from a 10x Chromium.
# This is the step in which you analyze differential expression between conditions, and then identify clusters that might be transcriptionally different.
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
#singleCellClusterMarkersPath <- args[4]
#demuxId <- args[5]

# Import data
#dataCombinedPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/from20190720_excitCluster/analyzeClusters/seuratObj_data.combined.RData'
load(dataCombinedPath)

#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
metadata<-read.csv(metadataPath)

#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/from20190720_excitCluster'

# singleCellClusterMarkersPath <- '/labs/tpalmer/resources/singleCellClusterMarkers/20190729_markerGenes_comprehensive.csv'
#singleCellClusterMarkers <- read.csv(singleCellClusterMarkersPath)

# mouse 16p region genes
mouse16pGenes <- c('Bola2', 'Sult1a1', 'Coro1a', 'Mapk3', 'Gdpd3', 'Ypel3', 'Tbx6', 'Ppp4c', 'Aldoa', 'Fam57b', 'Doc2a', 
                   'Ino80e', 'Hirip3', 'Taok2', 'Tmem219', 'Kctd13', 'Asphd1', 'Sez6l2', 'Cdipt', 'Mvp', '290009e17rik', 
                   'Prrt2', 'Maz', 'Kif22', '1810010M01Rik', 'Al467606', 'Qprt', 'Spn', 'Cd2bp2', 'Tbc1d10b', 'Mylpf',
                   'C16orf92', 'Bcl7c')

# Create folder to store output
setwd(outputDir)
subDir <- 'differentialExpression'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}


# # Import demultiplexed IDs if possible
# argsLen <- length(args);
# 
# outfile <- if (argsLen < 5) {
#   print('Using non-demultiplexed IDs...')
# }else{
#   demuxID <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/sexDemux/mapCellsToNewIDs.RData'
#   print(paste0('Using demultiplexed IDs from from ', demuxID))
#   load(demuxID)
#   
#   data.combined@meta.data$SampleID <- mapCellsToNewIDs[match(row.names(data.combined@meta.data), mapCellsToNewIDs$cells),'demuxID']
#   
#   # Generate some UMAPs on the basis of sex
#   bc_sexDemuxFolder <- 'batchEffects_sexDemux'
#   dir.create(file.path(outputDir, subDir, bc_sexDemuxFolder))
#   setwd(file.path(outputDir, subDir, bc_sexDemuxFolder))
#   
#   p1 <- DimPlot(data.combined, reduction = "umap", group.by = 'SampleID')
#   p_split <- DimPlot(data.combined, reduction = "umap", split.by = 'SampleID', ncol=3)
#     
#   ggsave(filename = paste0('umap_batch_demuxID.pdf'), plot = p1, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 20, height=20, units ='cm')
#   ggsave(filename = paste0('umap_batch_splitby_demuxID.pdf'), plot = p_split, device='pdf', path = file.path(outputDir, subDir, bc_sexDemuxFolder), width = 40, height=40, units ='cm')
#   
# }

data.combi <- data.combined

# Chose to run on assay=RNA (raw data):
# https://github.com/satijalab/seurat/issues/1198
# https://github.com/satijalab/seurat/issues/1168

findDEgenes <- function(data.combi, myCluster, factorToCompare, level_A, level_B){
  # Switch identity to the condition you want to compare across
  data.combi$comparison <- paste(data.combi@meta.data[,'seurat_clusters'], data.combi@meta.data[,factorToCompare], sep = "_")
  Idents(data.combi) <- "comparison"
  
  # Differential expresssion
  grp_A <- paste0(myCluster, '_', level_A)
  grp_B <- paste0(myCluster, '_', level_B)
  #deGenes_AvB <- FindMarkers(data.combi, ident.1 = grp_A, ident.2 = grp_B, verbose = FALSE, logfc.threshold = 0)
  #deGenes_AvB <- FindMarkers(data.combi, ident.1 = grp_A, ident.2 = grp_B, verbose = FALSE, logfc.threshold = 0, min.pct = 0)
  deGenes_AvB <- FindMarkers(data.combi, ident.1 = grp_A, ident.2 = grp_B, verbose = FALSE, logfc.threshold = 0, assay='RNA')
  head(deGenes_AvB, n = 15)
  
  deGenes_AvB_signif <- deGenes_AvB[which(deGenes_AvB$p_val_adj < 0.05),]
  
  write.table(deGenes_AvB, file = paste0('deGenes_clust', myCluster, '_', level_A, 'vs', level_B, '.csv'), quote=FALSE)
  write.table(deGenes_AvB_signif, file = paste0('deGenes_signif_clust', myCluster, '_', level_A, 'vs', level_B, '.csv'), quote=FALSE)
  
  return(deGenes_AvB_signif)
}


for (clust in levels(data.combined@meta.data$seurat_clusters)[c(1:3)]){
  print(paste0('Finding out differentially expressed genes for ', clust))
  
  # Compare between specific conditions
  print('Comparing specific conditions...')
  subDir_sc <- 'specificConditions'

  if (file.exists(file.path(subDir_sc))){
    setwd(file.path(outputDir, subDir, subDir_sc))
  } else {
    dir.create(file.path(outputDir, subDir, subDir_sc))
    setwd(file.path(outputDir, subDir, subDir_sc))
  }
  print(paste0('Setting path to ', getwd() ) )
  assign(paste0('de_WT.SALvWT.LPS_', clust), findDEgenes(data.combined, clust, 'cond', 'WT.SAL', 'WT.LPS') )
  assign(paste0('de_WT.SALvHET.SAL_', clust), findDEgenes(data.combined, clust, 'cond', 'WT.SAL', 'HET.SAL') )
  #assign(paste0('de_WT.SALvHET.LPS_', clust), findDEgenes(data.combined, clust, 'cond', 'WT.SAL', 'HET.LPS') )
  #assign(paste0('de_HET.SALvHET.LPS_', clust), findDEgenes(data.combined, clust, 'cond', 'HET.SAL', 'HET.LPS') )
  
  # # Compare along sex
  # print('Comparing along sex generally...')
  # subDir_sex <- 'Sex'
  # 
  # if (file.exists(file.path(subDir, subDir_sex))){
  #   setwd(file.path(outputDir, subDir, subDir_sex))
  # } else {
  #   dir.create(file.path(outputDir, subDir, subDir_sex))
  #   setwd(file.path(outputDir, subDir, subDir_sex))
  # 
  # }
  # 
  # print(paste0('Setting path to ', getwd() ) )
  # assign(paste0('de_sex_', clust), findDEgenes(data.combined, clust, 'Sex', 'M', 'F') )
  # 
  # # Compare along Genotype
  # print('Comparing along genotype generally...')
  # subDir_geno <- 'Genotype'
  # 
  # if (file.exists(file.path(subDir, subDir_geno))){
  #   setwd(file.path(outputDir, subDir, subDir_geno))
  # } else {
  #   dir.create(file.path(outputDir, subDir, subDir_geno))
  #   setwd(file.path(outputDir, subDir, subDir_geno))
  # 
  # }
  # 
  # print(paste0('Setting path to ', getwd() ) )
  # assign(paste0('de_geno_', clust), findDEgenes(data.combined, clust, 'Genotype', 'WT', 'HET') )
  # 
  # # Compare along Injection Condition
  # print('Comparing along injection condition generally...')
  # subDir_cond <- 'Condition'
  # 
  # if (file.exists(file.path(subDir, subDir_cond))){
  #   setwd(file.path(outputDir, subDir, subDir_cond))
  # } else {
  #   dir.create(file.path(outputDir, subDir, subDir_cond))
  #   setwd(file.path(outputDir, subDir, subDir_cond))
  # 
  # }
  # 
  # print(paste0('Setting path to ', getwd() ) )
  # assign(paste0('de_cond_', clust), findDEgenes(data.combined, clust, 'Condition', 'SAL', 'LPS') )

}

# Visualize features of interest

## 16p genes
FeaturePlot(data.combined, features = c("Ppp4c", "Aldoa", "Hirip3", 'Ypel3'), split.by = "cond", max.cutoff = 3,
            cols = c("grey", "red"))

FeaturePlot(data.combined, features = c("Coro1a", "Fam57b", 'Kif22'), split.by = "cond", max.cutoff = 3,
            cols = c("grey", "red"))

## WT.SAL vs. WT. LPS
FeaturePlot(data.combined, features = c("Eif3j1", 'Ubc', 'Hist1h2ap'), split.by = "cond", max.cutoff = 3,
            cols = c("grey", "red"))


# Bar plot - which clusters have the most DE genes? !!!!!!
filesToPlot_vWT.LPS <- ls(pattern='de_WT.SALvWT.LPS_*')
filesToPlot_vHET.SAL <- ls(pattern='de_WT.SALvHET.SAL_*')
filesToPlot_vHET.LPS <- ls(pattern='de_WT.SALvHET.LPS_*')

signifGenesCount <- data.frame(file=c(), nSignifGenes=c(), CtrlVersus=c())

# Count number of signif genes per cluster
for (f in filesToPlot_vWT.LPS){
  tmpRow <- data.frame(file = f, 
                       nSignifGenes = eval(parse(text=f)) %>% filter(p_val_adj < 0.05) %>% nrow,
                       ctrlVersus = 'WT.LPS',
                       clust = strsplit(f, split='_')[[1]][3])
  signifGenesCount <- rbind(signifGenesCount, tmpRow)
}

for (f in filesToPlot_vHET.SAL){
  tmpRow <- data.frame(file = f, 
                       nSignifGenes = eval(parse(text=f)) %>% filter(p_val_adj < 0.05) %>% nrow,
                       ctrlVersus = 'HET.SAL',
                       clust = strsplit(f, split='_')[[1]][3])
  signifGenesCount <- rbind(signifGenesCount, tmpRow)
}


for (f in filesToPlot_vHET.LPS){
  tmpRow <- data.frame(file = f, 
                       nSignifGenes = eval(parse(text=f)) %>% filter(p_val_adj < 0.05) %>% nrow,
                       ctrlVersus = 'HET.LPS',
                       clust = strsplit(f, split='_')[[1]][3])
  signifGenesCount <- rbind(signifGenesCount, tmpRow)
}


# Create plot
b <- ggplot(signifGenesCount, aes(x=clust, y=nSignifGenes, fill = ctrlVersus)) + 
  geom_bar(stat = 'identity') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0('bar_nSignifGenes_byClust.pdf'), plot = b, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 40, height=20, units ='cm')


# plot differences in means
plotDiffMeans <- function(data.combi, myCluster, identToUse, x_axis, y_axis, genesToLabel){

  # Subset out data of interest
  dataToPlot <- subset(data.combi, idents = myCluster)
  Idents(dataToPlot) <- identToUse
  avg.dataToPlot <- log1p(AverageExpression(dataToPlot, verbose = FALSE)$RNA)
  avg.dataToPlot$gene <- rownames(dataToPlot)

  # create plot
  p1 <- ggplot(avg.t.cells, aes_string(x_axis, y_axis)) + geom_point() +
    ggtitle(paste0("Scatterplot of average gene expression \nacross conditions ", x_axis, " and ", y_axis))
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  print(p1)
}

plotDiffMeans(data.combined, '9', 'cond', 'WT.SAL', 'WT.LPS', c('Xist') )




# Bar plot - what genes are DE in which proportion of clusters (according to type)

# Bar plot - what genes are EXPRESSED in which proportion of clusters (according to type)


# Violin plot visualization of 16p reion genes
plots <- VlnPlot(data.combined, features = mouse16pGenes, split.by = "stim", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)


#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
#save(ls(pattern='de_geno_*'), file = 'deGenes_byGenotype.RData')
save.image(file = paste0(paste0("allVars.RData")))
#save(data.combined, file = 'seuratObj_data.combined.RData')

print('~*~ All done! ~*~')









