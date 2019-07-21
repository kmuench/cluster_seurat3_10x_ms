# 1_initializeVars.R
# Kristin Muench
# 2019.07.16
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
#load(args[1])
mtxPath <- args[3]
metadataPath <- args[2]
outputDir <- args[1]
#altID_cells_path <- args[4]

# #Uncomment if you want to only focus on a subset of clusters
# whatClustUse <- c('4','2','0', '5', '12','6','14','10','3','9','1','22','20')

# !!! FOR TROUBLESHOOTING - UNCOMMENT FOR PATHS
#mtxPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/mtxData'
#metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190718_sampleTable.csv'
#outputDir <- '/scratch/users/kmuench/output/cnv16p/201907_cluster_seurat_10x_ms/20190716_remakePipeline/'
# altID_cells_path <- read.csv("/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel/allNewID_cells.csv")
# altID_cells <- read.csv("/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/originalSampleBarcodeLUT.csv")

# !!!!! START HERE

print('Run : 1_initializeVars.R')
print(paste0('Mtx File location: ', mtxPath))
print(paste0('Metadata location: ', metadataPath))
print(paste0('Output location: ', outputDir))
#print(paste0('Alternative Cell IDs location: ', altID_cells_path))

# # load cell name IDs (pooled or demultiplexed)
# altID_cells <- read.csv(altID_cells_path)

# make subdirectory
setwd(outputDir)
subDir <- 'initializeVars'

if (file.exists(subDir)){
  setwd(file.path(outputDir, subDir))
} else {
  dir.create(file.path(outputDir, subDir))
  setwd(file.path(outputDir, subDir))
}

# get file names
list.filenames <- list.files(mtxPath)

print('Example file location:')
print(paste0(mtxPath, '/', list.filenames[1],'/mm10/'))

# read metadata
metadata<-read.csv(metadataPath)

# create look up table for list.filenames with corresponding sample IDs?
files <- data.frame(filenames = list.filenames)
files$ids <- unlist(lapply(strsplit(list.filenames, "_"), '[[', 1))

# pull those sample IDs into a new list.filenames
list.filenames.wt.sal <- files[files$ids %in% metadata[metadata$Group == 'WT.SAL','SampleID'], 'filenames']
list.filenames.wt.lps <- files[files$ids %in% metadata[metadata$Group == 'WT.LPS','SampleID'], 'filenames']
list.filenames.het.sal <- files[files$ids %in% metadata[metadata$Group == 'HET.SAL','SampleID'], 'filenames']
list.filenames.het.lps <- files[files$ids %in% metadata[metadata$Group == 'HET.LPS','SampleID'], 'filenames']

print(files)
print(list.filenames.wt.sal)

# # load subclustering IDs
# load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/groupCompare/subclust_IDs_c4_2_0_5.RData') # interneurons and dienceph
# subclust_IDs <- tst


# Create Seurat object for each condition
## First sample imported to start the Seurat object
## Each subsequent loop is used to tack on Seurat objects corresponding to the other samples
## filenames: list of filenames with only files relevant to that condition, e.g. list.filenames.wt.sal
## mtxPath: path to the directory storing 10x mtx data
## condname: name of this condition, e.g. 'WT.SAL'
## files: Lookup table for filename and sampleID
myCondSeurat <- function(filenames, mtxPath, condName, files, minCells, nFeatures_chosen, nFeature_RNA_thresh){
  print('~*~')
  print(paste0('Working on ', condName))
  
  # pull out only files of interest
  files_relevant <- files[which(files$filenames %in% filenames),]
  print(files_relevant)
  
  # list of objects you'll make
  listOfObjects <- c()
  #allCellIDs <- c()
  
  # now merge with the others
  for (i in 1:length(filenames)){
    
    print(paste0('Making initial Seurat object with ', filenames[i], '...'))
    nameOfThisObject <- paste0('data_', filenames[i])
    print(paste0('nameOfThisObject = ', nameOfThisObject))
    listOfObjects <- c(listOfObjects, nameOfThisObject)
    
    print('Reading 10X data...')
    data_tmp_imported <- Read10X(data.dir = paste0(mtxPath, '/', filenames[i], '/mm10/') )
    print('Create Seurat Object...')
    data_tmp <- CreateSeuratObject(counts = data_tmp_imported, project = condName, min.cells = minCells)
    data_tmp$orig.ident <- files_relevant[i,'ids'] # Add sample IDs
    
    
    # filter out cells if desired for this round
    data_tmp <- subset(data_tmp, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500) # numbers need to be typed in
    
    
    # LIST OF CELLS TO REMOVE GOES HERE
    # subset(data_tmp, cell_id %in% (putative doublets))
    
    # normalize data and find variable featuers
    print('Normalize the data...')
    #data_tmp <- NormalizeData(data_tmp, verbose = FALSE)
    #secondData_obj <- NormalizeData(secondData_obj, verbose = FALSE)
    print('Find variable features data...')
    data_tmp <- FindVariableFeatures(data_tmp, selection.method = "vst", nfeatures = nFeatures_chosen)
    #secondData_obj <- FindVariableFeatures(secondData_obj, selection.method = "vst", nfeatures = nFeatures_chosen)
    
    print('Rename data_tmp...')
    assign(nameOfThisObject, data_tmp)
    
    #allCellIDs <- c(allCellIDs, files_relevant[i,2])
  }
  
  # merge this new seurat object with existing data
  print(paste0('Merging Seurat data...'))
  
  
  data <- merge(eval(parse(text=listOfObjects[1])), c(eval(parse(text=listOfObjects[2])), eval(parse(text=listOfObjects[3])), eval(parse(text=listOfObjects[4])) ), project = condName, 
                      add.cell.ids = files_relevant[,2], merge.data = TRUE )
  
  # give everything in this Seurat Object a condition name
  print('Giving everything in merged object a condition name...')
  data$cond <- condName
  
  # ## ** SUBCLUSTERING ** subset data according to subclustering - ONLY 13/25
  # data <- SubsetData(data, cells.use = row.names( subclust_IDs[ which(subclust_IDs$res.1.2 %in% whatClustUse) ,] ))
  
  # give everything an alternative ID
  print('Giving everything in merged object an alternative ID...')
  metadata_df <- data@meta.data
  data$trueBarcode <- sapply(strsplit(row.names(metadata_df), split="_"),tail, n=1L)
  #data$sample <- altID_cells[ match( row.names(data@meta.data), altID_cells$Barcode ) , 'sample']
  
  return(data)
}

# Do without any filtering of nFeature_RNA

wt.sal <- myCondSeurat(filenames = list.filenames.wt.sal, 
                       mtxPath = mtxPath, 
                       condName = 'WT.SAL', 
                       files = files, 
                       minCells = 3, 
                       nFeatures_chosen = 2000,
                       nFeature_RNA_thresh = TRUE)
wt.sal <- FindVariableFeatures(wt.sal, selection.method = "vst", nfeatures = 2000)

wt.lps <- myCondSeurat(filenames = list.filenames.wt.lps, 
                       mtxPath = mtxPath, 
                       condName = 'WT.LPS', 
                       files = files, 
                       minCells = 3, 
                       nFeatures_chosen = 2000,
                       nFeature_RNA_thresh = TRUE)
wt.lps <- FindVariableFeatures(wt.lps, selection.method = "vst", nfeatures = 2000)


het.sal <- myCondSeurat(filenames = list.filenames.het.sal, 
                       mtxPath = mtxPath, 
                       condName = 'HET.SAL', 
                       files = files, 
                       minCells = 3, 
                       nFeatures_chosen = 2000,
                       nFeature_RNA_thresh = TRUE)
het.sal <- FindVariableFeatures(het.sal, selection.method = "vst", nfeatures = 2000)


het.lps <- myCondSeurat(filenames = list.filenames.het.lps, 
                       mtxPath = mtxPath, 
                       condName = 'HET.LPS', 
                       files = files, 
                       minCells = 3, 
                       nFeatures_chosen = 2000,
                       nFeature_RNA_thresh = TRUE)
het.lps <- FindVariableFeatures(het.lps, selection.method = "vst", nfeatures = 2000)

# Print what is in each one
print('Samples included in each variable:')
table(wt.sal@meta.data$orig.ident)
table(wt.lps@meta.data$orig.ident)
table(het.sal@meta.data$orig.ident)
table(het.lps@meta.data$orig.ident)


# Perform integration
print('Performing integration...')
print('Finding integration anchors...')
data.anchors <- FindIntegrationAnchors(object.list = list(wt.sal, wt.lps, het.sal, het.lps), dims = 1:49)
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:49)

# Visualize Features for QC Purposes
## Run these after running myCondSeurat without visualizeFeatures(wt.sal)
## Visualize QC metrics as a violin plot to check what nFeature_RNA should be
## based on this, choose thresholding of 1000 (lower bound) and 7500 (upper bound)

## create subdirectory to store output
subDir_visFeatures <- 'visualizeFeatures'
setwd(file.path(outputDir, subDir))
if (file.exists(subDir_visFeatures)){
  setwd(file.path(outputDir, subDir,subDir_visFeatures))
} else {
  dir.create(file.path(outputDir, subDir, subDir_visFeatures))
  setwd(file.path(outputDir, subDir, subDir_visFeatures))
}

## generate feature vis graphs
visualizeFeatures <- function(myObject){
  v <- VlnPlot(myObject, features = c("nFeature_RNA", "nCount_RNA"))
  v1 <- VlnPlot(data.combined, features = c("nFeature_RNA", "nCount_RNA"), group.by='orig.ident')
  v2 <- VlnPlot(data.combined, features = c("nFeature_RNA", "nCount_RNA"), group.by='cond')
  f <- FeatureScatter(myObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  ggsave(filename = paste0(deparse(substitute(myObject)), '_vlnPlot.pdf'), plot = v, device='pdf', path = file.path(outputDir, subDir), width = 30, height=20, units ='cm')
  ggsave(filename = paste0(deparse(substitute(myObject)), '_vlnPlot_byOrigIdent.pdf'), plot = v1, device='pdf', path = file.path(outputDir, subDir), width = 30, height=20, units ='cm')
  ggsave(filename = paste0(deparse(substitute(myObject)), '_vlnPlot_byCond.pdf'), plot = v2, device='pdf', path = file.path(outputDir, subDir), width = 30, height=20, units ='cm')
  ggsave(filename = paste0(deparse(substitute(myObject)), '_featureScatter.pdf'), plot = f, device='pdf', path = file.path(outputDir, subDir), width = 30, height=20, units ='cm')
  
}

visualizeFeatures(wt.sal)
visualizeFeatures(wt.lps)
visualizeFeatures(het.sal)
visualizeFeatures(het.lps)
visualizeFeatures(data.combined)


# Use Jackstraw plots to determine the dimensionality of data and back-justify better dimensions for data anchors
## justification: estimated number using https://github.com/satijalab/seurat/issues/1248

## create subdirectory to store output
subDir_jackstraw <- 'jackstraw'
setwd(file.path(outputDir, subDir))
if (file.exists(subDir_jackstraw)){
  setwd(file.path(outputDir, subDir,subDir_jackstraw))
} else {
  dir.create(file.path(outputDir, subDir, subDir_jackstraw))
  setwd(file.path(outputDir, subDir, subDir_jackstraw))
}

## perform jackstraw for each individual object
for (ds in c('wt.sal', 'wt.lps', 'het.sal', 'het.lps')){
  ds_explore <- ScaleData(eval(parse(text=ds)), features = rownames(eval(parse(text=ds))) )
  ds_explore <- RunPCA(ds_explore, features = VariableFeatures(object = ds_explore), npcs=80) #!! ERROR: max(nu, nv) must be positive
  ds_explore <- JackStraw(ds_explore, num.replicate = 100, dims=80)
  ds_explore <- ScoreJackStraw(ds_explore, dims = 30:70) #looking for where sharp dropoff
  j <- JackStrawPlot(ds_explore, dims = 30:70)
  ggsave(filename = paste0(ds, '_jackstraw.tiff'), plot = j, device='tiff', path = file.path(outputDir, subDir), width = 30, height=12, units ='cm')
}


# save variables
print('saving variables...')
setwd(file.path(outputDir, subDir))
save.image(file = paste0("allVars.RData"))
save(wt.sal, file = 'seuratObj_wt.sal.RData')
save(wt.lps, file = 'seuratObj_wt.lps.RData')
save(het.sal, file = 'seuratObj_het.sal.RData')
save(het.lps, file = 'seuratObj_het.lps.RData')
save(data.combined, file = 'seuratObj_data.combined.RData')


print('~*~ All done! ~*~')







# ## BONUS ## MAKING NEW CSV FOR POOLED RUNTHRU
# 
# tst<- data.frame(Barcode = as.character(row.names(wt.sal@meta.data)), 
#                  sample = wt.sal@meta.data$orig.ident)
# tst[, ] <- lapply(tst[, ], as.character)
# tst2<- data.frame(Barcode = row.names(wt.lps@meta.data), sample = wt.lps@meta.data$orig.ident)
# tst2[, ] <- lapply(tst2[, ], as.character)
# tst3<- data.frame(Barcode = row.names(het.sal@meta.data), sample = het.sal@meta.data$orig.ident)
# tst3[, ] <- lapply(tst3[, ], as.character)
# tst4<- data.frame(Barcode = row.names(het.lps@meta.data), sample = het.lps@meta.data$orig.ident)
# tst4[, ] <- lapply(tst4[, ], as.character)
# 
# originalSampleBarcodeLUT <- rbind(tst,tst2,tst3,tst4)
# setwd('/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/')
# write.csv(originalSampleBarcodeLUT, file = 'originalSampleBarcodeLUT.csv', quote=FALSE, row.names=FALSE)
