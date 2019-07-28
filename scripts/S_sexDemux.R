# S_sexDemux.R
# Kristin Muench
# 2019.07.28
# Pipeline to run Seurat to analyze Mouse scRNA-Seq data from a 10x Chromium.
# This script analyzes expression of murine sex genes in dataset and demultiplexes them
#
#
# Explains the basic QC steps: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
#
# Following tutorial here: https://satijalab.org/seurat/v3.0/immune_alignment.html
# 
# Integration explained: https://satijalab.org/seurat/v3.0/integration.html
#
# Documentation of assay slot: https://github.com/satijalab/seurat/wiki/Assay
#
# # # # # # # # # # # # # # # # # #

# Load needed library
library('Seurat') 
library('cowplot')
library('ggplot2')

print(paste0('Seurat package Version: ', packageVersion('Seurat')) )
print(paste0('Cowplot package Version: ', packageVersion('cowplot')) )
print(paste0('Ggplot2 package Version: ', packageVersion('ggplot2')) )


# import data

## import paths from command line
args <- commandArgs(TRUE)
outputDir <- args[1]
metadataPath <- args[2]
dataCombinedPath <- args[3]
sexGenesPath <- args[4]

## for troubleshooting
# outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500'
# metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
# dataCombinedPath <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190720_s1s2_normalization_filter1000filter7500/cluster/allVars_res1.2.RData'
# sexGenesPath <- '/labs/tpalmer/projects/cnv16p/resources/sexGenes/mouse_sexGenes_20190728.csv'

## load metadata
metadata<-read.csv(metadataPath)
load(dataCombinedPath)
sexGenes <- read.csv(sexGenesPath)
sexGenes_uniqList <- as.character(unique(sexGenes$Gene.name))

## declare known samples
knownM <- c('KM2', 'KM13', 'KM1', 'J', 'KM18')  
knownF <- c('I', 'K', 'L')  
pooled <- c('H', 'M', 'F', 'N', 'E','O','G','P')

## declare known cells
cells_M <- row.names(data.combined@meta.data[data.combined@meta.data$orig.ident %in% knownM,])
cells_F <- row.names(data.combined@meta.data[data.combined@meta.data$orig.ident %in% knownF,])
cells_unk  <- row.names(data.combined@meta.data[data.combined@meta.data$orig.ident %in% pooled,])
  

# Create folder to store output
setwd(outputDir)
subDir <- 'sexDemux'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}



# Select discriminating genes to use

## Normalized data matrices
data_norm_knownSex <- subset(data.combined, 
                             cells = c(cells_M, cells_F) ) # note that I am plotting RAW COUNTS here
data_norm_knownSex@meta.data$Sex <- 'M'
data_norm_knownSex@meta.data[which(row.names(data_norm_knownSex@meta.data) %in% cells_F),'Sex'] <- 'F'

data_norm_M <- data_norm_knownSex[,colnames(data_norm_knownSex) %in% cells_M]
data_norm_F <- data_norm_knownSex[,colnames(data_norm_knownSex) %in% cells_F]
data_sexGenes <- data_norm_knownSex[which(row.names(data_norm_knownSex) %in% sexGenes_uniqList),]

## Find average gene expression for sex genes
geneAvgs <- data.frame(Male = rowMeans(data_norm_M), Female = rowMeans(data_norm_F))
geneAvgs$delta <- abs((geneAvgs$Male - geneAvgs$Female))/(geneAvgs$Male + geneAvgs$Female) # which genes have the largest scaled difference between the M and F average?
sexGeneAvgs <- geneAvgs[which(row.names(geneAvgs) %in% sexGenes_uniqList),] # only look at sex genes

## visualize gene candidates in the subset of data that is for known samples
### Older gene candidates - uncomment to view
# plots <- VlnPlot(data_norm_knownSex, features = c("Uty", 'H2afb1','Magix', 'Ddx3y','Tsix','Kdm5d'), split.by = "Sex", group.by = "cond",
#         pt.size = 0, combine = FALSE, slot = 'RNA')
# CombinePlots(plots = plots, ncol = 1)
# 
# plots2 <- VlnPlot(data_norm_knownSex, features = c("Xist", "Eif2s3y", "Eif2s3x", "Hmgb3", "Tmsb4x"), split.by = "Sex", group.by = "orig.ident",
#         pt.size = 0, combine = FALSE, slot = 'RNA')
# CombinePlots(plots = plots2, ncol = 1)

### Choose genes that look most promising
chosenSexGenes <- c('Xist', 'Ddx3y', 'Eif2s3y')

demoPlots <- VlnPlot(data_norm_knownSex, features = chosenSexGenes, split.by = "Sex", group.by = "orig.ident",
                     pt.size = 0, combine = FALSE, assay = 'RNA') # include assay flag to use non-normalized data
demoPlots_toSave <- CombinePlots(plots = demoPlots, ncol = 3)
ggsave(filename = paste0(paste0('vlnPlot_chosenSexGenes.pdf')), 
       plot = demoPlots_toSave, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 45, height=15, units ='cm')

## visualize expression of X-linked genes in tSNE - X = genes, Y=sample
data.combined_small <- subset(data.combined, 
                              cells = row.names(data.combined@meta.data[data.combined@meta.data$orig.ident %in% c('KM2','J','I','K','H','M'),]) ) # note that I am plotting RAW COUNTS here
fp_chosenSexGenes <- FeaturePlot(data.combined_small, features = chosenSexGenes, split.by = "orig.ident")
ggsave(filename = paste0(paste0('fp_chosenSexGenes.pdf')), 
       plot = fp_chosenSexGenes, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 20, height=20, units ='cm')

## ROC curve to select expression threshold for Chosen Genes


# Characterize those sex-linked genes

## In singlets, how many cells are...

### Xist+, Y1+ AND Y2+?
### Xist+, Y1+ OR Y2+?
### Xist+, Y1- AND Y2-?
### Xist-, Y1+ OR Y2+?
### Xist-, Y1+ OR Y2+?
### Xist-, Y1- AND Y2-?
### Xist-, Y1- OR Y2-?


# Determine what typical expression is in singlets

## Violin plot - Xist vs other two

## What number of cells in each range?

# Create factor for whether or not something is expressed




## Of the F cells that are Xist-, 

### what categories do they fall in?

### how often do they express Y marker genes?

# Compare that to what is expressed in doublets

# Remove cells expressing X and Y

## Where are they in UMAP? (entire cluster need to be removed?)

## List of cells to remove

# Generate list of M and F samples

# Rename samples, redo datafame?