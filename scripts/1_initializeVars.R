# initializeVars.R
# Kristin Muench
# 2019.07.16
# Pipeline to run Seurat setup for mouse scRNA-Seq data with CCA
# # # # # # # # # # # # # # # # # #

# Load needed library
# Load older version of Seurat which has the functions I need
library('Seurat', lib.loc='/scg/apps/software/r/3.5.0/scg/seurat_2.3')
packageVersion('Seurat')

## run this in case you accidentally load Seurat 3
#detach("package:Seurat", unload=TRUE)


# import from command line
args <- commandArgs(TRUE)
#load(args[1])
mtxPath <- args[1]
metadataPath <- args[2]
outputDir <- args[3]
altID_cells_path <- args[4]
whatClustUse <- c('4','2','0', '5', '12','6','14','10','3','9','1','22','20')

# # !!! FOR TROUBLESHOOTING - UNCOMMENT FOR PATHS
# mtxPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/mtxData'
# metadataPath <- '/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/20190211_sampleTable.csv'
# outputDir <- '/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190228_runThroughDemultiplex/'
#altID_cells_path <- read.csv("/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190219_makeSparse/20190226_sexLabel/allNewID_cells.csv")
#altID_cells <- read.csv("/labs/tpalmer/projects/cnv16p/data/scRNASeq/mouse/metadata/originalSampleBarcodeLUT.csv")

print('Run CCA Track: cca_makeVars.R')
print(paste0('Mtx File location: ', mtxPath))
print(paste0('Metadata location: ', metadataPath))
print(paste0('Output location: ', outputDir))
print(paste0('Alternative Cell IDs location: ', altID_cells_path))

# load cell name IDs (pooled or demultiplexed)
altID_cells <- read.csv(altID_cells_path)

# make subdirectory
setwd(outputDir)
subDir <- 'makeVars'

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

# # load subclustering IDs
# load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/groupCompare/subclust_IDs_c4_2_0_5.RData') # interneurons and dienceph
# subclust_IDs <- tst

# for when you are fed up and just want to import everything!
load('/scratch/users/kmuench/output/cnv16p/201901_cluster_pooled_10x_ms/20190329_troubleshootRunthrough_demux/clusters/data.combined_withClust_r1.2_CC40.RData')
subclust_IDs <- data.combined@meta.data

# Create Seurat object for each condition
## First sample imported to start the Seurat object
## Each subsequent loop is used to tack on Seurat objects corresponding to the other samples
## filenames: list of filenames with only files relevant to that condition, e.g. list.filenames.wt.sal
## mtxPath: path to the directory storing 10x mtx data
## condname: name of this condition, e.g. 'WT.SAL'
## files: Lookup table for filename and sampleID
myCondSeurat <- function(filenames, mtxPath, condName, files){
  print('~*~')
  print(paste0('Working on ', condName))
  
  # pull out only files of interest
  files_relevant <- files[which(files$filenames %in% filenames),]
  print(files_relevant)
  
  # make first seurat object
  print(paste0('Making initial Seurat object with ', filenames[1], '...'))
  initialData <- Read10X(data.dir = paste0(mtxPath, '/', filenames[1], '/mm10/') )
  data <- CreateSeuratObject(raw.data = initialData, project = condName, min.cells = 3, min.features=200, names.field = files_relevant[1,2])
  
  ## Add sample IDs
  data@meta.data$orig.ident <- files_relevant[1,'ids'] 
  
  # merge with the second
  # make this iteration's seurat object
  print(paste0('Reading in data from ',mtxPath, '/', filenames[2],'...'))
  secondData <- Read10X(data.dir = paste0(mtxPath, '/', filenames[2], '/mm10/') )
  secondData_obj <- CreateSeuratObject(raw.data = secondData, project = condName, min.cells = 3, min.features=200, names.field = files_relevant[2,2]) # create object
  secondData_obj@meta.data$orig.ident <- files_relevant[2,'ids']
  
  # merge this new seurat object with existing data
  print(paste0('Merging Seurat data with data from ', filenames[2], '...'))
  data <- MergeSeurat(data, secondData_obj, project = condName, 
                      add.cell.id1 = files_relevant[1,2], 
                      add.cell.id2 = files_relevant[2,2])
  
  # now merge with the others
  for (i in 3:length(filenames)){
    
    # make this iteration's seurat object
    print(paste0('Reading in data from ',mtxPath, '/', filenames[i],'...'))
    addingData <- Read10X(data.dir = paste0(mtxPath, '/', filenames[i], '/mm10/') )
    addingData_obj <- CreateSeuratObject(raw.data = addingData, project = condName, min.cells = 3, min.features=200, names.field = files_relevant[i,2]) # create object
    addingData_obj@meta.data$orig.ident <- files_relevant[i,'ids']
    
    # merge this new seurat object with existing data
    print(paste0('Merging Seurat data with data from ', filenames[i], '...'))
    data <- MergeSeurat(data,addingData_obj, project = condName, add.cell.id2 = files_relevant[i,2])
  }
  
  # give everything in this Seurat Object a condition name
  data@meta.data$cond <- condName
  
  ## ** SUBCLUSTERING ** subset data according to subclustering - ONLY 13/25
  data <- SubsetData(data, cells.use = row.names( subclust_IDs[ which(subclust_IDs$res.1.2 %in% whatClustUse) ,] ))
  
  # prefilter now that all the samples are in
  data <- FilterCells(data, subset.names = "nGene", low.thresholds = 200, high.thresholds = Inf)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
  data <- ScaleData(data, display.progress = F, vars.to.regress = "nUMI", do.par = TRUE, num.cores = 4) # took out genes.use = hv.genes - seems complicated with separate seurat obj
  
  # give everything an alternative ID
  data@meta.data$trueBarcode <- sapply(strsplit(row.names(data@meta.data ), split="_"),tail, n=1L)
  data@meta.data$sample <- altID_cells[ match( row.names(data@meta.data), altID_cells$Barcode ) , 'sample']
  
  return(data)
}


wt.sal <- myCondSeurat(list.filenames.wt.sal, mtxPath, 'WT.SAL', files)
wt.lps <- myCondSeurat(list.filenames.wt.lps, mtxPath, 'WT.LPS', files)
het.sal <- myCondSeurat(list.filenames.het.sal, mtxPath, 'HET.SAL', files)
het.lps <- myCondSeurat(list.filenames.het.lps, mtxPath, 'HET.LPS', files)

# Print what is in each one
table(wt.sal@meta.data$orig.ident)
table(wt.lps@meta.data$orig.ident)
table(het.sal@meta.data$orig.ident)
table(het.lps@meta.data$orig.ident)

# save variables
setwd(file.path(outputDir, subDir))
save.image(file = paste0("allVars.RData"))
save(wt.sal, file = 'seuratObj_wt.sal.RData')
save(wt.lps, file = 'seuratObj_wt.lps.RData')
save(het.sal, file = 'seuratObj_het.sal.RData')
save(het.lps, file = 'seuratObj_het.lps.RData')

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
