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
library('dplyr')

print(paste0('Seurat package Version: ', packageVersion('Seurat')) )
print(paste0('Cowplot package Version: ', packageVersion('cowplot')) )
print(paste0('Ggplot2 package Version: ', packageVersion('ggplot2')) )


# import data
print('Beginning R Script, importing data...')


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
print('Selecting discriminating genes to use...')

## Normalized data matrices
data_norm_knownSex <- subset(data.combined, 
                             cells = c(cells_M, cells_F) ) # note that I am plotting RAW COUNTS here
data_norm_knownSex@meta.data$Sex <- 'M'
data_norm_knownSex@meta.data[which(row.names(data_norm_knownSex@meta.data) %in% cells_F),'Sex'] <- 'F'

data_norm_knownSex_exprs <- FetchData(data_norm_knownSex, vars = sexGenes_uniqList)
data_norm_M <- data_norm_knownSex_exprs[rownames(data_norm_knownSex_exprs) %in% cells_M,]
data_norm_F <- data_norm_knownSex_exprs[rownames(data_norm_knownSex_exprs) %in% cells_F,]
#data_sexGenes <- data_norm_knownSex[which(row.names(data_norm_knownSex) %in% sexGenes_uniqList),]

## Find average gene expression for sex genes
geneAvgs <- data.frame(Male = colMeans(data_norm_M), Female = colMeans(data_norm_F))
geneAvgs$delta <- abs((geneAvgs$Male - geneAvgs$Female))/(geneAvgs$Male + geneAvgs$Female) # which genes have the largest scaled difference between the M and F average?
#sexGeneAvgs <- geneAvgs[which(row.names(geneAvgs) %in% sexGenes_uniqList),] # only look at sex genes

## visualize gene candidates in the subset of data that is for known samples
### Older gene candidates - uncomment to view
# plots <- VlnPlot(data_norm_knownSex, features = c("Uty", 'H2afb1','Magix', 'Ddx3y','Tsix','Kdm5d'), split.by = "Sex", group.by = "cond",
#         pt.size = 0, combine = FALSE, slot = 'RNA')
# CombinePlots(plots = plots, ncol = 1)
# 
# plots2 <- VlnPlot(data_norm_knownSex, features = c("Xist", "Eif2s3y", "Eif2s3x", "Hmgb3", "Tmsb4x"), split.by = "Sex", group.by = "orig.ident",
#         pt.size = 0, combine = FALSE, slot = 'RNA')
# CombinePlots(plots = plots2, ncol = 1)

## Choose genes that look most promising
chosenSexGenes <- c('Xist', 'Ddx3y', 'Eif2s3y')

## Visualize in violin plot what those discriminating genes look like
demoPlots <- VlnPlot(data_norm_knownSex, features = chosenSexGenes, split.by = "Sex", group.by = "orig.ident",
                     pt.size = 0, combine = FALSE, assay = 'RNA') # include assay flag to use non-normalized data
demoPlots_toSave <- CombinePlots(plots = demoPlots, ncol = 3)
ggsave(filename = paste0(paste0('vlnPlot_chosenSexGenes.pdf')), 
       plot = demoPlots_toSave, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 45, height=15, units ='cm')

# !!! NOT WORKING
# ## visualize expression of X-linked genes in tSNE - X = genes, Y=sample
# data.combined_small <- subset(data.combined, 
#                               cells = row.names(data.combined@meta.data[data.combined@meta.data$orig.ident %in% c('I','K','KM13','KM18','M'),]) ) # note that I am plotting RAW COUNTS here
# fp_chosenSexGenes <- FeaturePlot(data.combined[['RNA']], features = chosenSexGenes, split.by = "orig.ident", slot = "counts")
# 
# ggsave(filename = paste0(paste0('fp_chosenSexGenes.pdf')), 
#        plot = fp_chosenSexGenes, device='pdf', 
#        path = file.path(outputDir, subDir), 
#        width = 40, height=20, units ='cm')


## ROC curve to select expression threshold for Chosen Genes 
print('Plotting ROC curve...')
library('plotROC')

### prepare ROC curve plot data
plotROC <- data.frame(cellName = row.names(data_norm_knownSex@meta.data), Sex = data_norm_knownSex@meta.data$Sex)
chosenSexGene_exprs <- FetchData(data_norm_knownSex[['RNA']], vars = chosenSexGenes, slot = "counts")
plotROC <- merge(plotROC, chosenSexGene_exprs, by.x = 'cellName', by.y = 'row.names')
plotROC$binarySexLabel <- 0
plotROC[which(plotROC$Sex == 'F'), 'binarySexLabel'] <- 1

longroc <- melt_roc(plotROC, "binarySexLabel", c("Xist", "Eif2s3y", "Ddx3y"))

plotThreeROC <- ggplot(longroc, aes(d = D, m = M, color = name)) + geom_roc(cutoffs.at = c(0, 0.5, 1, 3, 7, 13)) + style_roc()
ggsave(filename = paste0(paste0('ROC_curve.pdf')), 
       plot = plotThreeROC, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 25, height=15, units ='cm')

### choose a threshold
thresh_Xist <- 1
thresh_Eif2s3y <- 0
thresh_Ddx3y <-0


# Characterize cells with low and high Xist in the known sex (single sample) cells
print('Characterizing the cells with low and high Xist in the singly-sequenced samples...')

plotROC$y1_or_y2 <- 'FALSE'
plotROC[which(plotROC$Ddx3y > thresh_Ddx3y | plotROC$Eif2s3y > thresh_Eif2s3y),'y1_or_y2'] <- 'TRUE'

plotROC$y1_and_y2 <- 'FALSE'
plotROC[which( (plotROC$Ddx3y > thresh_Ddx3y) & (plotROC$Eif2s3y > thresh_Eif2s3y) ),'y1_and_y2'] <- 'TRUE'

## How many cells have each combination of markers?
### NOTE: WORKS BEST WHEN THE Y-MARKERS ARE <= AND XIST IS >
plotROC %>% group_by(Sex) %>% count(Xist>thresh_Xist)
plotROC %>% group_by(Sex) %>% count(Ddx3y>thresh_Ddx3y)
plotROC %>% group_by(Sex) %>% count(Eif2s3y>thresh_Eif2s3y)

## Count how many cells fall into each possible combination of Xist, Eif2s3y, and Ddx3y expression
count_bins <- plotROC %>% group_by(Sex) %>% count(Xist>thresh_Xist,y1_and_y2, y1_or_y2)
write.csv(count_bins, file = paste0('counts_knownSexSamples_xistThresh', thresh_Xist, '_EifThresh', thresh_Eif2s3y, '_DdxThresh', thresh_Ddx3y, '.csv') )


## Of the cells that are triple negative, what cell types do they fall under?
### We want to know if we call all triple negative cells Male which female cell types that will affect
tripleNegCells <- data.frame(cellName = plotROC[ which( !((plotROC$Xist > thresh_Xist)) & (plotROC$y1_or_y2 == 'FALSE') & (plotROC$Sex == 'F') ) , 'cellName'] )
tripleNegCells <- data.combined@meta.data[match( tripleNegCells$cellName, row.names(data.combined@meta.data) ) , c("seurat_clusters", "CellsPerSample", "nFeature_RNA") ]

tripleNegCells <- tripleNegCells %>% mutate(nFeature_RNA_category=cut(nFeature_RNA, breaks=c(-Inf, 2045, 3053, Inf), labels=c("low","middle","high")))


### make pie chart
#### tips: http://www.sthda.com/english/wiki/ggplot2-pie-chart-quick-start-guide-r-software-and-data-visualization
makePieChart <- function(tripleNegCells, clusterBy){
  # prep data
  pieData <- tripleNegCells %>% group_by_at(clusterBy) %>% tally()
  pieData$percent <- (pieData$n / sum(pieData$n) )*100
  
  pieData$labels <- apply(pieData[,clusterBy], 2, function(x) paste0(x, ' [', round(pieData$percent, 2), '%]') )
  
  pieData <- pieData[rev(order(pieData$percent)),]
  pieData$labels <- factor(pieData$labels, levels = pieData$labels)
  
  # blank theme for pie chart
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  # make pie
  p <- ggplot(pieData, aes_string(x=factor(1), y='percent', fill='labels')) +
    geom_bar(width = 1, stat = "identity", color = "black") + 
    coord_polar("y", start=0) +
    theme_minimal() +
    theme(axis.text.x=element_blank()) +
    blank_theme
    
    # add when you have colors decided:
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  
  print(p)
  ggsave(filename = paste0(paste0('pieChart_F_tripleNeg_groupby', clusterBy, '.pdf')), 
         plot = p, device='pdf', 
         path = file.path(outputDir, subDir), 
         width = 30, height=15, units ='cm')
}

makePieChart(tripleNegCells, 'seurat_clusters')
makePieChart(tripleNegCells, 'CellsPerSample')
makePieChart(tripleNegCells, 'nFeature_RNA_category') # so most of these cells in the lowest third of expression



# Apply this scheme to the pooled samples
print('Applying scheme to all pooled samples...')

## cell IDs with Xist > 1
exprs <- FetchData(data.combined[['RNA']], vars = chosenSexGenes, slot = "counts")
cells_above_thresh_Xist <-row.names(exprs[which(exprs$Xist > thresh_Xist),])

data.combined@meta.data$Sex <- 'M'
data.combined@meta.data[which(row.names(data.combined@meta.data) %in% cells_above_thresh_Xist),'Sex'] <- 'F'
data.combined@meta.data[which(row.names(data.combined@meta.data) %in% cells_F), 'Sex'] <- 'F' # for those 102 ambiguous cells, can take advantage of the fact that we already know label

data.combined@meta.data$SampleID <- apply(data.combined@meta.data, 1, function(x) paste0(x['SampleID'], '_', x['Sex']))

## Are there any cells that have both Xist above thresh and at least one Y-markers above thresh?
### As we would expect, only F cells have these doubles (because everything with high Xist called F), and only pooled cells have detected doubles
exprs$y1_or_y2 <- 'FALSE'
exprs[which(exprs$Ddx3y > thresh_Ddx3y | exprs$Eif2s3y > thresh_Eif2s3y),'y1_or_y2'] <- 'TRUE'

exprs$XandY <- 'FALSE'
exprs[which( (exprs$y1_or_y2 == TRUE) & exprs$Xist > thresh_Xist),'XandY'] <- 'TRUE'

cells_XandY <- row.names(exprs[which(exprs$XandY == 'TRUE'),])

data.combined@meta.data$XandY <- 'FALSE'
data.combined@meta.data[row.names(data.combined@meta.data) %in% cells_XandY, 'XandY'] <- 'TRUE'

## Visualize where those cells appear
XandY_loc_all <- DimPlot(data.combined, cells.highlight = cells_XandY)
XandY_loc_bySampleID <- DimPlot(data.combined, cells.highlight = cells_XandY, split.by = 'SampleID', ncol=4)

ggsave(filename = paste0(paste0('dp_cellsExpressingXandY.pdf')), 
       plot = XandY_loc_all, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 30, height=30, units ='cm')

ggsave(filename = paste0(paste0('dp_cellsExpressingXandY_bySampleID.pdf')), 
       plot = XandY_loc_bySampleID, device='pdf', 
       path = file.path(outputDir, subDir), 
       width = 30, height=30, units ='cm')

## !!!!! ADD IN BAR PLOT PROPORTION VISUALIZATION

#  Save variables
print('Saving variables...')
setwd(file.path(outputDir, subDir))
save.image(file = paste0(paste0("allVars.RData")))
save(data.combined, file = 'seuratObj_data.combined_demux.RData')


print('~*~ All done! ~*~')

