###########    -----   Get data to start for the    ------    ###########
###########    -- DICE Tissue - Healthy Tissue Project  --    ###########

# Version.
# Full version: 1.0

# Whole version: 1.0
# Copy from version 0.7 where we have triplechecked the content and prepared the script to be reported on our repo.


### -------------------------- Description -------------------------- ###
# Script to get the basic data to get the final figures for the DICE Tissue project.

# NOTE: If any of the objects previously processed need any further processing, you may as well load the already processed object and do so for this session.


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(stringr)
library(data.table)
library(tidyr)
library(gtools)


### --------------------------- Functions --------------------------- ###
coll.version <- 0.4
gen.data.path <- '/path/to/this/dir' # Path where all preprocessed files have been saved. These will be the input to produce the final figures.
coll.file <- paste0(gen.data.path, '/functions_collection.', coll.version, '.R') # This file is also provided on the repo (same dir)
file.exists(coll.file)
source(coll.file)


### ----------------------- General Arguments ----------------------- ###

# ---> General definitions.
# @ Seed.
set.seed(seed=1)
# @ Cluster labels for each dataset.
obj.extended.names <- c(
  hlty.cd4='Healthy-CD4',
  hlty.cd8='Healthy-CD8',
  hlty.nk='Healthy-NK',
  hlty.b='Healthy-B',
  hlty.mye='Healthy-Myeloid'
)
clust.labs <- c(
  hlty.cd4='RNA_snn_res.0.3',
  hlty.cd8='RNA_snn_res.0.2',
  hlty.nk='RNA_snn_res.0.1',
  hlty.b='RNA_snn_res.0.1',
  hlty.mye='RNA_snn_res.0.2'
)
# ---> Path definitions.
# Module signatures.
module.sigs.path <- paste0(gen.data.path, '/module_signatures')
# Output directory to save the list of cell barcodes kept acroos separate cell type-specific analyses.
custom.subset.dir <- '/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/custom_subsets/all_cell-types_aggr'
# Individual aggr table file paths
aggr.table.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-14-2022/aggr/data'
aggr.data.path <- '/mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr'
# ---> File definitions.
# Seurat objects.
obj.files <- c(
  hlty.cd4='/path/to/results/from/processing/CD4/CD4_Round-2/output/seurat_objects/SeuratObjectForPrjCD4_Round-2_WithArgs_NoPCs_30.RDS',
  hlty.cd8='/path/to/results/from/processing/CD8/CD8_Round-2/output/seurat_objects/SeuratObjectForPrjCD8_Round-2_WithArgs_NoPCs_30.RDS',
  hlty.nk='/path/to/results/from/processing/NK/NK_Round-5/output/seurat_objects/SeuratObjectForPrjNK_Round-5_WithArgs_NoPCs_30.RDS',
  hlty.b='/path/to/results/from/processing/B/B_Round-4/output/seurat_objects/SeuratObjectForPrjB_Round-4_WithArgs_NoPCs_30.RDS',
  hlty.mye='/path/to/results/from/processing/Myeloid/Myeloid_Round-3/output/seurat_objects/SeuratObjectForPrjMyeloid_Round-3_WithArgs_NoPCs_30.RDS'
)
# Demuxlet results assessment.
demuxlet.reports.path <- '/path/to/results/from/demuxlet_assessment' # FILES TO BE PROVIDED. AVAILABLE UPON REQUEST? @KEVIN TO DETERMINE IF THEY'RE NEEDED DURING PEER REVIEW
# Donor metadata file.
donor.meta.file <- '/path/to/this/dir/DonorMetadataFromCristian.0.1.csv'
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(module.sigs.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', module.sigs.path))
reports.path <- gen.data.path
essential.files <- obj.files
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# ---> Seurat obejcts.
obj.list <- lapply(X=obj.files, FUN=readRDS)
# ---> Donor metadata file.
donor.meta <- fread(file=donor.meta.file)

# ---> Define aggr tables.
# Define cell and tissue types.
aggr.pffx.lab <- 'R24_Cancer_Batches-1-to-20_'
cell.types <- c('CD4', 'CD8', 'Myeloid', 'NKAndB')
tissue.types <- c('Normal-lung', 'Tumor-lung')
aggr.files <- paste0(aggr.table.path, '/', aggr.pffx.lab, cell.types)
aggr.files <- paste0(rep(x=aggr.files, each=2), '_', tissue.types)
aggr.files <- paste0(aggr.files, '_aggr_table_annotated.1.0.csv')
if(!all(file.exists(aggr.files))) stop('Unexpected error.\n')
names(aggr.files) <- paste0(aggr.pffx.lab, rep(x=cell.types, each=2), '_', tissue.types)


### ---------------------- Data preprocessing ----------------------- ###

# ---> Donor Metadata
donor.meta[, donor.id.tag:=str_replace(string=ID, pattern='DLCP[0]*', replacement='')]
donor.meta[, ID:=NULL]
donor.meta[, age.tag:=Age]; donor.meta[, Age:=NULL]
donor.meta[, sex.tag:=Gender]; donor.meta[, Gender:=NULL]
donor.meta[, performance.status.tag:=`Performance Status`]; donor.meta[, `Performance Status`:=NULL]
donor.meta[, smoking.status.tag:=`Smoking Status`]; donor.meta[, `Smoking Status`:=NULL]
donor.meta[, smoking.info.tag:=`Smoking info`]; donor.meta[, `Smoking info`:=NULL]
donor.meta[, ext.smoking.status.tag:=`Smoking Status Extended`]; donor.meta[, `Smoking Status Extended`:=NULL]

# ---> Confirm whether seurat objects have been previously preprocessed. If so, avoid repeating the process again.
tmp.files <- paste0(reports.path, '/SeuratObj_', obj.extended.names, '.RDS')
names(tmp.files) <- names(obj.extended.names)
# For now, keep only the healthy tissue data. Note that this step must be erased eventually.
tmp.files <- grep(x=tmp.files, pattern='Healthy', value=TRUE)
tmp.check <- all(file.exists(file=tmp.files))


if(!tmp.check){
  # ---> Apply preprocessing.
  obj.list <- lapply(X=names(obj.list), FUN=function(tmp.obj) preprocess.obj(seurat.obj=obj.list[[tmp.obj]], set.lab=tmp.obj, clusts.lab=clust.labs[tmp.obj], donor.meta.data=donor.meta, get.metadata=FALSE))
  names(obj.list) <- names(obj.files)

  # ---> Specific removal of cell clusters.
  # ---> HLTY CD4
  tmp.lab <- 'hlty.cd4'
  # Clusters to remove: 8
  tmp.remove <- '8'
  tmp.data <- obj.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
  tmp.cells <- Cells(obj.list[[tmp.lab]])[!tmp.data %in% tmp.remove]
  obj.list[[tmp.lab]] <- subset(x=obj.list[[tmp.lab]], cells=tmp.cells)
  # ---> HLTY CD8
  tmp.lab <- 'hlty.cd8'
  # Clusters to remove: 7 & 8.
  tmp.remove <- c('7', '8', '9')
  tmp.data <- obj.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
  tmp.cells <- Cells(obj.list[[tmp.lab]])[!tmp.data %in% tmp.remove]
  obj.list[[tmp.lab]] <- subset(x=obj.list[[tmp.lab]], cells=tmp.cells)
  # ---> HLTY NK
  tmp.lab <- 'hlty.nk'
  # Clusters to remove: None.
  # ---> HLTY B
  tmp.lab <- 'hlty.b'
  # Clusters to remove: 4.
  tmp.remove <- c('4')
  tmp.data <- obj.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
  tmp.cells <- Cells(obj.list[[tmp.lab]])[!tmp.data %in% tmp.remove]
  obj.list[[tmp.lab]] <- subset(x=obj.list[[tmp.lab]], cells=tmp.cells)
  # ---> HLTY Myeloid.
  tmp.lab <- 'hlty.mye'
  # Clusters to remove: 8.
  tmp.remove <- c('8')
  tmp.data <- obj.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
  tmp.cells <- Cells(obj.list[[tmp.lab]])[!tmp.data %in% tmp.remove]
  obj.list[[tmp.lab]] <- subset(x=obj.list[[tmp.lab]], cells=tmp.cells)

  # ---> Save seurat objects.
  lapply(X=names(obj.list), FUN=function(tmp.obj){
    tmp.file.name <- paste0(reports.path, '/SeuratObj_', obj.extended.names[tmp.obj], '.RDS')
    if(file.exists(tmp.file.name)) saveRDS(file=tmp.file.name, object=obj.list[[tmp.obj]]) else cat(paste0('File already exists for object: ', tmp.obj, '\n'))
  })
}else{
  obj.list <- lapply(X=tmp.files, FUN=readRDS)
  names(obj.list) <- names(tmp.files)
}