############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary figure stack    -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# Version.
# Full version: 2.0

# Whole version: 2.0
# Copy from version 1.7 where we have triplechecked the content and prepared the script to be reported on our repo.


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Script to get the final figures for the DICE Tissue project.
# Some files used in here might not be provided in repo, but might be available upon request. See contact info in the general repo or at the DICE website.


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
library(data.table)
library(tidyr)
library(Seurat)
library(stringr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(magrittr)


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############
coll.version <- 0.4
gen.data.path <- '/path/to/this/dir' # Path where all preprocessed files have been saved. These will be the input to produce the final figures.
coll.file <- paste0(gen.data.path, '/functions_collection.', coll.version, '.R') # This file is also provided on the repo (same dir)
file.exists(coll.file)
source(coll.file)


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# @ Paper.
main.prj <- 'R24'
this.prj <- 'R24_Cancer'
this.figure <- 'figure_stack_v1'
tissue.type <- 'healthy'
# @ Seed.
set.seed(seed=1)
# @ Subset size.
subset.n <- 50000
# @ Dataset labels.
obj.extended.names <- c(
  hlty.cd4='Healthy-CD4',
  hlty.cd8='Healthy-CD8',
  hlty.nk='Healthy-NK',
  hlty.b='Healthy-B',
  hlty.mye='Healthy-Myeloid'
)
# @ Cluster labels for each dataset.
clust.labs <- c(
  hlty.cd4='RNA_snn_res.0.3',
  hlty.cd8='RNA_snn_res.0.2',
  hlty.nk='RNA_snn_res.0.1',
  hlty.b='RNA_snn_res.0.1',
  hlty.mye='RNA_snn_res.0.2'
)
# @ For plotting purposes.
blank.complement.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_blank(), axis.line=element_blank()) # Default blank.
blank.complement.2 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=0.4), axis.line=element_line(size=0.4), axis.ticks.length=unit(0.15, "cm")) # For the dot plot.
# @ Gene signature definitions.
cd4.ctl.sign <- 'cell.cytotoxicity.patil.score'
# @ For tag analysis.
# min.no.cells.per.group <- 50 # @ Internal parameter for tag specific analysis (part of figure ?).
colors.table <- NULL
# ---> Cluster identities
clusts.defs <- list(
  `hlty.cd4`=c(
    `0`='TRM',
    `1`='TRM/TCM',
    `2`='TCM',
    `3`='Conv-CTLs',
    `4`='TREG',
    `5`='NonConv-CTLs',
    `6`='TFH',
    `7`='THIFNR'
  ),
  `hlty.cd8`=c(
    `0`='TEFF',
    `1`='TRM',
    `2`='GZMKhigh',
    `3`='TCM',
    `4`='KIR+',
    `5`='MAIT',
    `6`='THIFNR'
  ),
  `hlty.nk`=c(
    `0`='CD56dim/CD16+',
    `1`='CD56bright/CD16neg'
  ),
  `hlty.b`=c(
    `0`='0',
    `1`='1',
    `2`='Plasmablasts',
    `3`='BRM'
  ),
  `hlty.mye`=c(
    `0`='CM',
    `1`='NCM',
    `2`='TR-Macrophages',
    `3`='Macrophage/cDC2',
    `4`='cDC2',
    `5`='pDC',
    `6`='cDC1',
    `7`='mregDC',
    `8`='MAST'
  )
)
# ---> Color defintions.
# @ Color per cell type.
cell.type.cols <- c(
  Myeloid='#B22222',
  NK='#FFD700',
  B='#32CD32',
  CD4='#00BFFF',
  CD8='#EE82EE'
)
# @ Color per cluster per dataset.
# See definitions at: https://docs.google.com/spreadsheets/d/1_lNI4M7P-LBuwfQWEgag4pHaycMZMNcmzGv_k_e7x5U/edit#gid=0
clusts.cols <- list(
  `hlty.cd4`=c(
    `0`='#3333FF', # TRMs
    `1`='#00BFFF', # TRM w/ circulating features.
    `2`='#BFEFFF', # TCM
    `3`='#E60000', # Conventional CD4-CTLs
    `4`='#319272', # TREG
    `5`='#FF9933', # Unconventional CD4-CTLs
    `6`='#A6B0E3', # TFH
    `7`='#FFCC33' # THIFNR
  ),
  `hlty.cd8`=c(
    `0`='#CC99FF', # TEFF
    `1`='#8000FF', # TRM
    `2`='#4DA6FF', # GZMK+
    `3`='#FF1AB3', # TCM
    `4`='#39AC39', # Innate activating/inhibiting receptor
    `5`='#676765', # MAIT-like
    `6`='#FFCC33' # THIFNR
  ),
  `hlty.nk`=c(
    `0`='#FFD700', # CD56dim/CD16+
    `1`='#E87600' # CD56bright/CD16-
  ),
  `hlty.b`=c(
    `0`='#32CD32', # Undefined
    `1`='#196719',
    `2`='#3366FF', # Plasmablasts
    `3`='#E6B3FF'
  ),
  `hlty.mye`=c(
    `0`='#B22222', # Classical monocytes
    `1`='#FF8C00', # Non-classical monocytes
    `2`='#0052CC', # Tissue-resident macrophages
    `3`='#80BFFF', # Macrophage/cDC2 signature mixture
    `4`='#EFA9A9', # cDC2
    `5`='#BFBFBF', # pDC
    `6`='#803E96', # cDC1
    `7`='#9EC026', # mregDC
    `8`='#A6A6A6' # MAST
  )
)
# @ Other terms' colors.
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# ---> Path definitions.
prj.gen.path <- '/path/to/this/dir' # USER TO CHANGE
gen.data.path <- paste0(prj.gen.path, '/paper_items')
final.figs.path <- paste0(prj.gen.path, '/final_figures')
this.fig.path <- paste0(final.figs.path, '/', this.figure)
if(!dir.exists(this.fig.path)) dir.create(this.fig.path)
reports.path <- paste0(this.fig.path, '/', this.figure, '_panels')
# ---> File definitions.
# Seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', obj.extended.names, '.RDS')
names(objs.files.list) <- names(obj.extended.names)
# Features info.
features.info.file <- paste0(gen.data.path, '/FeaturesInfo.tsv')
# Specific files.
#     eGene results for individual datasets.
sg.egene.res.file <- paste0(gen.data.path, '/eGeneResults_Hlty.csv') # Only available upon request.
#     eGene results for whole sets.
wh.egene.res.file <- paste0(gen.data.path, '/eGeneResults-WholeSets_Hlty.csv') # Only available upon request.
#     MSigDB signatures.
msigdb.file <- paste0(gen.data.path, '/msigdb/MSigDB_HALLMARK_2023-03-31.csv')
#     Donors' metadata.
# donor.meta.file <- paste0(gen.data.path, '/DonorMetadata.csv')
#     Sex-biased genes.
sex.bias.genes.file <- paste0(gen.data.path, '/SexBiasedTranscriptDiscovery_Strategy-PerClusterPerCellType.csv')
#     Donors' cell counts from FACS data.
facs.counts.file <- paste0(gen.data.path, '/FACSCellCountInfo.csv')
#     Raw and normalized cell counts per donor.
cell.counts.file <- paste0(gen.data.path, '/ExtendedCellCountsPerDonor.RDS')
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(this.fig.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', this.fig.path))
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  lapply(X=objs.files.list, FUN=function(x) return(x))
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))
# ---> Present arguments (to do)


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------- Data Loading -------------------------- ###
############    -----------------------------------------    ############

# ---> Seurat obejcts.
# Main seurat objects.
srt.objs.list <- lapply(X=objs.files.list, FUN=readRDS)
names(srt.objs.list) <- names(objs.files.list)

# ---> Specific files.
#     Reference.
feature.info <- read.delim(file=features.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
#     eGene results.
sg.egene.res <- fread(file=sg.egene.res.file)
wh.egene.res <- fread(file=wh.egene.res.file)
#     MSigDB signatures.
msigdb <- fread(file=msigdb.file)
#     Donors' metadata.
# donor.meta.data <- fread(file=donor.meta.file)
#     Sex-biased genes.
sex.bias.genes <- fread(file=sex.bias.genes.file)
# ---> FACS cells counts.
facs.counts <- fread(file=facs.counts.file)
#     Raw and normalized cell counts per donor.
cell.pop.counts <- readRDS(file=cell.counts.file)


############    --- --------------------------------------    ############
### ---------------------- Data preprocessing ----------------------- ###
############    -----------------------------------------    ############

# ---> HLTY CD4
tmp.lab <- 'hlty.cd4'
# Clusters to combine: None.

# ---> HLTY CD8
tmp.lab <- 'hlty.cd8'
# Clusters to combine: None.

# ---> HLTY NK
tmp.lab <- 'hlty.nk'
# Clusters to combine: 0 & 2.
tmp.comb <- c('0', '2')
tmp.new.val <- tmp.comb[1]
tmp.data <- srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
tmp.data[tmp.data %in% tmp.comb] <- tmp.new.val
srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- tmp.data

# ---> HLTY B
tmp.lab <- 'hlty.b'
# Clusters to combine: None.

# ---> HLTY Myeloid.
tmp.lab <- 'hlty.mye'
# Clusters to combine: None.

# ---> Specific files.


# ---> Obtain cell subset of each dataset that's equal across cell types.
subset.list <- lapply(X=srt.objs.list, FUN=function(seurat.obj){
    sample(x=Cells(seurat.obj), size=subset.n)
})


############    -----------------------------------------    ############
### --------------------------- Figure 1 ---------------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.1.path <- paste0(reports.path, '/figure_1')
if(!dir.exists(fig.1.path)) dir.create(fig.1.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4


### ------------------------- Main Figure 1 ------------------------- ###

# ---> UMAP plots depicting clusters.
# Whole dataset.
for(tmp.obj in names(srt.objs.list)){
    meta.data <- as.data.table(cbind(
    srt.objs.list[[tmp.obj]]@meta.data,
    srt.objs.list[[tmp.obj]]@reductions$umap@cell.embeddings
    ))
    tmp.ggplot <- ggplot(data=meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clust.labs[tmp.obj])) +
    geom_point(alpha=1, size=0.5) +
    scale_color_manual(values=clusts.cols[[tmp.obj]]) +
    labs(x='UMAP 1', y='UMAP 2', col='', fill='Cluster')
    # Output.
    tmp.lab <- paste0(obj.extended.names[tmp.obj], '_PopsOnUMAP_Whole')
    publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name=tmp.lab, type='tiff', blank.comp=blank.complement.1, do.legend=FALSE)
}
# Subset
for(tmp.obj in names(srt.objs.list)){
    cell.subset <- subset.list[[tmp.obj]]
    meta.data <- as.data.table(cbind(
    srt.objs.list[[tmp.obj]]@meta.data[cell.subset, ],
    srt.objs.list[[tmp.obj]]@reductions$umap@cell.embeddings[cell.subset, ]
    ))
    tmp.ggplot <- ggplot(data=meta.data, aes_string(x='UMAP_1', y='UMAP_2', col=clust.labs[tmp.obj])) +
    geom_point(alpha=1, size=2) +
    scale_color_manual(values=clusts.cols[[tmp.obj]]) +
    labs(x='UMAP 1', y='UMAP 2', col='', fill='Cluster')
    # Output.
    tmp.lab <- paste0(obj.extended.names[tmp.obj], '_PopsOnUMAP_Subset')
    publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name=tmp.lab, type='tiff', blank.comp=blank.complement.1, do.legend=FALSE)
}
# Number of cells per immune population.
tmp.data <- sapply(X=names(srt.objs.list), FUN=function(tmp.obj){
  length(Cells(srt.objs.list[[tmp.obj]]))
})
tmp.data <- as.data.table(tmp.data); colnames(tmp.data) <- 'cell.count'
tmp.data[, cell.pop:=names(srt.objs.list)]
tmp.file.name <- paste0(fig.1.path, '/CellCounts.csv')
fwrite(x=tmp.data, file=tmp.file.name)


### --------------------- Supplementary Figures --------------------- ###
### --------------------- Associated w/ Figure 1 -------------------- ###

# ---> Dot plots.
# Define genes of interest.
these.markers <- list(
  `hlty.cd4`=c(
    'TCF7', 'CCR7', 'SELL', # TCM
    # 'GPR183', 'IL7R', 'IFNGR1', 'KLRB1', # TRM/TCM
    'HOPX', 'IFNG', 'CCL4', # TRM
    'GNLY', 'PRF1', 'GZMB', # Classic CD4-CTL
    'GZMK', 'GZMA', 'GZMM', # Non-classic CD4-CTL
    'FOXP3', 'TIGIT', # TREG
    'CTLA4', 'PDCD1', 'CXCL13', # # TFH
    'IFI6', 'MX1' # THIFNR
  ),
  `hlty.cd8`=c(
    'GZMK', 'CRTAM', # GZMK-expressing
    'PRF1', 'GNLY', 'FCGR3A', # Effector
    'KLRC2', 'XCL1', 'KIR2DL4', # KIR-expressing
    'ZNF683', 'LGALS3', 'ITGAE', # TRM
    'TCF7', 'SELL', 'S1PR1', # TCM
    'KLRB1', 'TRAV1-2', 'TRBV6-4', # MAIT
    'IFI6', 'MX1' # THIFNR
  ),
  `hlty.b`=c(
    'CD82', 'TNFRSF13B',
    'FOXP1', 'YBX3',
    'ZEB2', 'FCRL3', 'FCRL5', 'CD19', 'CD72',
    'MZB1', 'IGHA1', 'IGHA2'
  ),
  `hlty.nk`=c(
    'FCGR3A', 'NKG7', 'GNLY', 'PRF1', 'GZMB', 'CX3CR1', # CD16+ NK cells.
    'XCL1', 'XCL2', 'ZNF683', 'GZMK' # CD16- NK cells.
  ),
  `hlty.mye`=c(
    'S100A9', 'S100A8', 'CD14', # CM
    'FCGR3A', 'CDKN1C', # NCM
    'OLR1', 'MRC1', 'MARCO', 'C1QC', # Macrophages
    'FCER1A', 'CD1C', 'IL1R2',
    'WDFY4', 'IRF8', 'CLEC9A',
    'IRF7', 'GZMB', 'TCF4',
    'LAMP3', 'BIRC3',  'CCR7'
    )
)
to.check.1 <- unlist(unique(these.markers))
to.check.1 <- translate.ids(ids=to.check.1, ensembl=FALSE)
to.check.1[is.na(names(to.check.1))]
# Define clusters' order.
these.pops <- list(
  `hlty.cd4`=c('2', '1', '0', '3', '5', '4', '6', '7'),
  `hlty.cd8`=c('2', '0', '4', '1', '3', '5', '6'),
  `hlty.b`=c('0', '1', '3', '2'),
  `hlty.nk`=c('0', '1'),
  `hlty.mye`=c('0', '1', '2', '3', '4', '6', '5', '7', '8')
)
# File heights.
file.heights <- c(
  `hlty.cd4`=8,
  `hlty.cd8`=8,
  `hlty.b`=5,
  `hlty.nk`=6,
  `hlty.mye`=8
)
# Plot per cell type.
# for(tmp.obj in names(srt.objs.list)){
for(tmp.obj in 'hlty.nk'){
  # if(tmp.obj=='hlty.b') next
  # Translate.
  tmp.markers <- translate.ids(ids=these.markers[[tmp.obj]], ensembl=FALSE)
  tmp.markers <- names(tmp.markers)
  names(tmp.markers) <- these.markers[[tmp.obj]]
  # Get plot.
  tmp.ggplot <- dot.plot(
    seurat.obj=srt.objs.list[[tmp.obj]], features=tmp.markers, slot='data', do.norm=FALSE, ensembl=TRUE,
    groups.tag=clust.labs[tmp.obj], groups.order=these.pops[[tmp.obj]], groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
    this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
    scale.by='radius', dot.scale=12, size.min=NA, size.max=NA,
    file.name=NULL
  )
  tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
  # Output.
  tmp.width <- (length(these.pops[[tmp.obj]]) * 0.5) + 1
  # tmp.height <- file.heights[tmp.obj]
  tmp.height <- (length(tmp.markers) * 0.4) + 1
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.1.path, file.name=paste0(obj.extended.names[tmp.obj], '_Markers_Opt-A'), type='pdf', blank.comp=blank.complement.2, do.legend=TRUE, width=tmp.width, height=tmp.height)
}


############    -----------------------------------------    ############
### --------------------------- Figure 2 ---------------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.2.path <- paste0(reports.path, '/figure_2')
if(!dir.exists(fig.2.path)) dir.create(fig.2.path)


### ------------------------- Main Figure 2 ------------------------- ###

# ---> Heatmap for IFN signature
# Do process for several gene sets.
# sig.names <- c('INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE', 'FATTY_ACID_METABOLISM', 'DNA_REPAIR', 'HYPOXIA')
sig.names <- c('INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE')
egene.labs <- c('Whole', 'Individual')
for(sig.name in sig.names){
  for(egene.lab in egene.labs){
    # Retrieve FDRs for genes in set.
    tmp.data.1 <- msigdb[name==sig.name & ensembl_in_ref, .(sig.name=name, ensembl=ensembl_gene, gene.name=gene_symbol)]
    tmp.data.2 <- if(egene.lab=='Whole') wh.egene.res[, .(ensembl=Geneid, cell.type, cluster, FDR, is.eGene)] else sg.egene.res[, .(ensembl=Geneid, cell.type, cluster, FDR, is.eGene)]
    tmp.data.1[!ensembl %in% tmp.data.2[, ensembl], uniqueN(ensembl)]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='ensembl', all.x=FALSE, all.y=FALSE)
    tmp.data[, id:=paste(cell.type, cluster, sep='.')]
    tmp.data <- tmp.data[, .(ensembl, id, FDR, is.eGene)]
    # Show overlap.
    data.to.plot <- list(
      tmp.data.1[ensembl %in% wh.egene.res[, Geneid], unique(ensembl)],
      tmp.data.2[is.eGene==TRUE, unique(ensembl)]
    )
    names(data.to.plot) <- c(str_replace_all(string=sig.name, pattern='_', replacement='\n'), 'eGene Set')
    whole.pop <- unique(intersect(x=msigdb[, ensembl_gene], wh.egene.res[, Geneid]))
    tmp.file.name <- paste0(fig.2.path, '/Ovlp_Set-', egene.lab, '_ForSig-', sig.name, '.tiff')
    venn.diagram(
      x=data.to.plot,
      hyper.test=TRUE,
      total.population=length(whole.pop),
      filename=tmp.file.name
    )
    # Retrieve data to plot.
    mat.to.plot <- copy(x=tmp.data)
    mat.to.plot[is.eGene==FALSE, FDR:=NA]; mat.to.plot[, is.eGene:=NULL]
    mat.to.plot <- spread(data=mat.to.plot, key=id, value=FDR, fill=NA)
    # Clean rows and columns.
    mat.to.plot <- as.data.frame(mat.to.plot, stringsAsFactors=FALSE)
    row.names(mat.to.plot) <- mat.to.plot$ensembl; mat.to.plot$ensembl <- NULL
    mat.to.plot <- as.matrix(mat.to.plot)
    genes.to.keep <- rowSums(mat.to.plot, na.rm=TRUE)>0
    pops.to.keep <- colSums(mat.to.plot, na.rm=TRUE)>0
    mat.to.plot <- mat.to.plot[genes.to.keep, pops.to.keep]
    # Retrieve data to apply clustering of the genes.
    tmp.data[is.eGene==TRUE, FDR:=1] # NOTE: COMMENT/UNCOMMENT THIS TO TAKE INTO CONSIDERATION FDRs FOR NON-eGENE.
    tmp.data[, is.eGene:=NULL]
    mat.to.clust <- spread(data=tmp.data, key=id, value=FDR, fill=NA)
    mat.to.clust <- as.data.frame(mat.to.clust, stringsAsFactors=FALSE)
    row.names(mat.to.clust) <- mat.to.clust$ensembl; mat.to.clust$ensembl <- NULL
    mat.to.clust <- as.matrix(mat.to.clust)
    tmp.check <- any(is.na(mat.to.clust)); if(tmp.check) stop('Unexpected error.\n')
    mat.to.clust <- mat.to.clust[genes.to.keep, pops.to.keep]
    # Apply clustering.
    mat.to.clust <- scale(x=mat.to.clust)
    dist.mat <- dist(x=mat.to.clust, method='euclidean')
    hclust.obj <- hclust(d=dist.mat, method='average')
    gene.order <- hclust.obj$labels[hclust.obj$order]
    # tmp.file.name <- paste0(fig.2.path, '/HClustProof.pdf')
    # pdf(file=tmp.file.name)
    # plot(hclust.obj)
    # dev.off()
    # Heatmap metadata.
    mat.to.plot <- mat.to.plot[gene.order, ]
    mat.to.plot <- -log10(mat.to.plot)
    mat.to.plot[mat.to.plot>10] <- 10
    tmp.meta <- data.frame(
        row.names=colnames(mat.to.plot),
        cell.type=str_extract(string=colnames(mat.to.plot), pattern='^[^\\.]+')
    )
    tmp.cols <- list(
      cell.type=cell.type.cols
    )
    # Adjust gene names.
    tmp.data.1 <- data.frame(mat.order=row.names(mat.to.plot))
    tmp.data.2 <- feature.info[, c('ensembl', 'name')]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by.x='mat.order', by.y='ensembl', sort=FALSE, all.x=TRUE, all.y=FALSE)
    if(any(is.na(tmp.data$name)) | nrow(tmp.data)!=nrow(mat.to.plot)) stop('Unexpected error.\n')
    row.names(mat.to.plot) <- tmp.data$name
    # Output.
    tmp.file.name <- paste0(fig.2.path, '/FDRHeatmap_Set-', egene.lab, '_ForSig-', sig.name, '.pdf')
    # pheatmap(mat=mat.to.plot, scale='row', cluster_rows=TRUE, cluster_cols=FALSE, filename=tmp.file.name)
    pheatmap(mat=mat.to.plot, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=tmp.meta, annotation_colors=tmp.cols, na_col='#808080', filename=tmp.file.name)
  }
}
files.to.rm <- list.files(path=fig.2.path, pattern='\\.log$', full.names=TRUE)
file.remove(files.to.rm)


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### --------------------- sex-biased transcripts -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.sbt.path <- paste0(reports.path, '/figure_on_sex_bias')
if(!dir.exists(fig.sbt.path)) dir.create(fig.sbt.path)


### -------------------------- Main Figure -------------------------- ###

# ---> Examples.

sbt.eg.path <- paste0(fig.sbt.path, '/to_choose')
if(!dir.exists(sbt.eg.path)) dir.create(sbt.eg.path)

tmp.data <- lapply(X=names(srt.objs.list), FUN=function(data.set){
  meta.data <- as.data.table(srt.objs.list[[data.set]]@meta.data)
  meta.data[, cluster:=get(clust.labs[data.set])]
  uniq.clusts <- meta.data[, sort(as.character(unique(cluster)))]
  # Define cell type-specific criteria
  cell.type <- str_replace(string=data.set, pattern='tumor\\.|hlty\\.', replacement='')
  tmp.report <- lapply(X=uniq.clusts, FUN=function(tmp.cluster){
    # Determine donors that have less than 10 cells for a given cluster.
    tmp.data <- meta.data[cluster==tmp.cluster, .N, by=full.donor.id.tag][N<=10, full.donor.id.tag]
    # Report on donor count for donors with < 10 cells.
    to.return <- data.table(
      donor.count=length(tmp.data),
      donors=paste0(tmp.data, collapse=';'),
      cluster.size=meta.data[cluster==tmp.cluster, .N]
    )
    return(to.return)
  })
  names(tmp.report) <- uniq.clusts
  tmp.report <- rbindlist(l=tmp.report, use.names=TRUE, idcol='cluster')
  return(tmp.report)
})
names(tmp.data) <- names(srt.objs.list)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='data.set')
tmp.data[donors=='', donors:=NA]


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### --------------------- phenotype associations -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.phen.assoc.path <- paste0(reports.path, '/figure_on_phen_associations')
if(!dir.exists(fig.phen.assoc.path)) dir.create(fig.phen.assoc.path)


### -------------------------- Main Figure -------------------------- ###

# ---> Function to retrieve population fractions

retrieve.immune.phens <- function(do.cd45.based=FALSE, do.logit=FALSE, report.dn=FALSE){
  # ---> Retrieve data.
    tmp.props <- rbindlist(l=cell.pop.counts, use.names=TRUE, idcol='data.set')
    # Add cluster description
    tmp.data <- lapply(X=clusts.defs, FUN=function(x){
        tmp.data <- data.table(
        pop.tag=names(x),
        clust.def=x
        )
        return(tmp.data)
    })
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='data.set')
    tmp.props <- merge(x=tmp.props, y=tmp.data, by=c('data.set', 'pop.tag'))
    tmp.props[,
        immune.subset:=paste(
            toupper(
                str_replace(string=data.set, pattern=ifelse(test=tissue.type=='healthy', yes='hlty\\.', no='tumor\\.'), replacement='')
            ),
            pop.tag,
            clust.def
        )
    ]
    # Type of fractions.
    if(do.cd45.based){
        tmp.props <- tmp.props[,
            .(donor, immune.subset, subset.value=subset.pseudo.frac)
        ]
        if(report.dn){
          tmp.props.2 <- facs.counts[tissue.type==ifelse(test=tissue.type=='healthy', yes='CL', no='CT'), .(donor=as.character(donor.id.tag), `DN T cells`=dn.cell.percent/100)]
          tmp.props <- merge(x=tmp.props, y=tmp.props.2, by='donor', all.x=TRUE, all.y=FALSE)
        }
    }else{
        # Sanity check.
        # tmp.data.1 <- tmp.props[, sum(subset.count), by=.(donor, data.set)]
        # tmp.data.2 <- unique(tmp.props[, .(donor, data.set, insil.count)])
        # tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor', 'data.set'))
        # tmp.data[V1!=insil.count]
        # That shows that I can go ahead and use the counts under column "insil.count"
        tmp.props <- tmp.props[,
            .(donor, immune.subset, subset.value=subset.count/insil.count)
        ] 
    }
    # Normalization
    if(do.logit) tmp.props[, subset.value:=car::logit(p=subset.value, percents=FALSE)]
    # Final wrangling
    tmp.props <- spread(data=tmp.props, key=immune.subset, value=subset.value, fill=0)
    # Return
    return(tmp.props)
}

# ---> Function to obtain the correlation plots.

among.pops.associations <- function(do.cd45.based=FALSE, do.logit=FALSE, sex.to.keep, report.dn=FALSE, tmp.reports.path){
    # @ Retrieve immune phenotypes.
    tmp.props <- retrieve.immune.phens(do.cd45.based=do.cd45.based, do.logit=do.logit, report.dn=report.dn)
    # @ Retrieve sex info and subset accordingly.
    tmp.data <- as.data.table(srt.objs.list[[1]]@meta.data)
    tmp.data <- unique(tmp.data[!is.na(full.donor.id.tag), .(donor=full.donor.id.tag, sex.tag)])
    tmp.props <- merge(x=tmp.props, y=tmp.data, by='donor')
    tmp.props <- tmp.props[sex.tag %in% sex.to.keep]
    tmp.props[, sex.tag:=NULL]
    # @ Sort columns according to cell type.
    col.order <- c(
      'donor',
      grep(x=colnames(tmp.props), pattern='^MYE', value=TRUE),
      grep(x=colnames(tmp.props), pattern='^NK', value=TRUE),
      grep(x=colnames(tmp.props), pattern='^B', value=TRUE),
      grep(x=colnames(tmp.props), pattern='^CD4', value=TRUE),
      grep(x=colnames(tmp.props), pattern='^CD8', value=TRUE),
      grep(x=colnames(tmp.props), pattern='^DN', value=TRUE)
    )
    tmp.check <- setdiff(x=col.order, y=colnames(tmp.props))
    if(length(tmp.check) > 0) stop('Unexpected error.\n')
    tmp.props <- tmp.props[, ..col.order]
    # @ Perform correlation analysis.
    corr.data <- rcorr(as.matrix(tmp.props)[, 2:ncol(tmp.props)], type='spearman')
    # Plot.
    tmp.file.name <- paste0(
        fig.phen.assoc.path,
        '/CorrAnalysis_Sex-', paste0(sex.to.keep, collapse='-'),
        '_CD45Normalized-', do.cd45.based,
        '_LogitTransformed-', do.logit, '.C.pdf'
    )
    tmp.r <- corr.data$r
    diag(tmp.r) <- NA
    tmp.p <- corr.data$P
    diag(tmp.p) <- NA
    pdf(file=tmp.file.name)
      corrplot(
          tmp.r,
          # type='upper',
          col=colorRampPalette(c("blue","white","red"))(200),
          order="original", diag=TRUE,
          p.mat=tmp.p, outline=FALSE,
          sig.level=0.05, insig='blank',
          addgrid.col='black', tl.col='black', shade.col='black', na.label.col='black'
      )
    dev.off()
    tmp.file.name <- paste0(
        fig.phen.assoc.path,
        '/CorrAnalysis_Sex-', paste0(sex.to.keep, collapse='-'),
        '_CD45Normalized-', do.cd45.based,
        '_LogitTransformed-', do.logit, '.B.pdf'
    )
    pdf(file=tmp.file.name)
      corrplot(
          corr.data$r,
          col=colorRampPalette(c("blue","white","red"))(200),
          order="original", diag=FALSE,
          p.mat=corr.data$P,
          sig.level=0.05, insig='blank',
          tl.pos='n', cl.pos='n'
      )
    dev.off()
    # Save reports.
    tmp.file.name <- paste0(
        fig.phen.assoc.path,
        '/CorrAnalysis_PVals',
        '_Sex-', paste0(sex.to.keep, collapse='-'),
        '_CD45Normalized-', do.cd45.based,
        '_LogitTransformed-', do.logit,
        '.csv'
    )
    write.csv(file=tmp.file.name, x=corr.data$P)
    tmp.file.name <- paste0(
        fig.phen.assoc.path,
        '/CorrAnalysis_CorrCoefs',
        '_Sex-', paste0(sex.to.keep, collapse='-'),
        '_CD45Normalized-', do.cd45.based,
        '_LogitTransformed-', do.logit,
        '.csv'
    )
    write.csv(file=tmp.file.name, x=corr.data$r)
    # @ Perform heirarchical clustering for cell subsets based on immune phenotypes
    # data.to.cluster <- corr.data$r
    # dist.mat <- as.dist(sqrt(2*(1-data.to.cluster))) # As suggested here: https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
    # # dist.mat <- dist(x=data.to.cluster, method='euclidean')
    # hclust.avg <- hclust(d=dist.mat, method='average')
    # # cut.avg <- cutree(tree=hclust.avg, k=6)
    # dend.obj <- as.dendrogram(hclust.avg)
    # dend.obj <- color_branches(dend=dend.obj, k=10, groupLabels=TRUE)
    # tmp.file.name <- paste0(
    #     fig.phen.assoc.path,
    #     '/HClust_Method-Avg',
    #     '_Sex-', paste0(sex.to.keep, collapse='-'),
    #     '_CD45Normalized-', do.cd45.based,
    #     '_LogitTransformed-', do.logit,
    #     '.pdf'
    # )
    # pdf(file=tmp.file.name, height=12, width=14)
    # plot(
    #   dend.obj, main=paste0('Pseudofractions-', do.cd45.based, ' Logit-', do.logit)
    # )
    # dev.off()
    # ---> Return correlation results.
    return(corr.data)
}

# Perform all possible types of analysis.
male.corr.res <- among.pops.associations(do.cd45.based=TRUE, do.logit=TRUE, sex.to.keep='Male', tmp.reports.path=fig.phen.assoc.path)
female.corr.res <- among.pops.associations(do.cd45.based=TRUE, do.logit=TRUE, sex.to.keep='Female', tmp.reports.path=fig.phen.assoc.path)
among.pops.associations(do.cd45.based=TRUE, do.logit=TRUE, sex.to.keep=c('Male', 'Female'), tmp.reports.path=fig.phen.assoc.path)

# ---> Comparison for the results between males and females.
# @ Prepare comparison matrix.
tmp.data.1 <- female.corr.res$r
diag(tmp.data.1) <- 100
tmp.data.1[female.corr.res$P > 0.05] <- NA
tmp.data.1 <- gdata::unmatrix(tmp.data.1)
# Code below was used as confirmation of how casting works between vectors and matrices.
# tmp.data <- as.vector(tmp.data.1)
# tmp.data.2 <- matrix(data=tmp.data, nrow=nrow(tmp.data.1), byrow=TRUE)
# all(tmp.data.1==tmp.data.2, na.rm=TRUE)
tmp.data.2 <- male.corr.res$r
diag(tmp.data.2) <- 100
tmp.data.2[male.corr.res$P > 0.05] <- NA
tmp.data.2 <- gdata::unmatrix(tmp.data.2)
# Condense
tmp.data <- data.table(
  association=names(tmp.data.1),
  to.check=names(tmp.data.2),
  female.r=tmp.data.1,
  male.r=tmp.data.2
)
# Sanity check.
tmp.check <- tmp.data[, all(association==to.check)]
if(tmp.check){
  tmp.data[, to.check:=NULL]
}else{
  stop('Unexpected error.\n')
}
# Proper format.
tmp.data <- separate(data=tmp.data, col=association, into=c('cell.type.1', 'cell.type.2'), sep=':')
tmp.data$cell.type.1 <- factor(x=as.character(tmp.data$cell.type.1), levels=colnames(male.corr.res$r))
tmp.data$cell.type.2 <- factor(x=as.character(tmp.data$cell.type.2), levels=rev(colnames(male.corr.res$r)))
# @ Define clases.
tmp.data[, class:='Unknown'] # Assume class isn't known.
#     Consistently non-significant associations 
tmp.data[is.na(female.r) & is.na(male.r), class:='Non-significant']
#     Consistently positive associations
tmp.data[female.r>0 & male.r>0, class:='Consistent, positive']
#     Consistently negative associations
tmp.data[female.r<0 & male.r<0, class:='Consistent, negative']
#     Inconsistent, only female
tmp.data[!is.na(female.r) & is.na(male.r), class:='Only female']
#     Inconsistent, only male
tmp.data[is.na(female.r) & !is.na(male.r), class:='Only male']
#     Diagonal values.
tmp.data[female.r==100 & male.r==100, class:='Diagonal']
#     At this point, all fields in the matrix have been covered. However, there could have also been inconsistent classes, where an association is significant for both sexes, yet positive for one and negative for the other.
# @ Define classes colors
tmp.cols <- c(
  `Non-significant`='#808080',
  `Consistent, positive`='#ff0000',
  `Consistent, negative`='#0000ff',
  `Only female`='#9370db',
  `Only male`='#e6e600',
  `Diagonal`='#ffffff'
)
# @ Report as heatmap.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cell.type.1, y=cell.type.2, fill=class)) +
  geom_tile() +
  scale_fill_manual(values=tmp.cols) +
  labs(x='', y='', fill='Class')
tmp.lab <- 'CorrAnalysisComp'
publish.plot(tmp.ggplot=tmp.ggplot, output.path=fig.phen.assoc.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=FALSE, do.rotate=TRUE)


############    -----------------------------------------    ############
### --------------------------- Figure X ---------------------------- ###
### --------------------- Further feature plots --------------------- ###
### ------------------------ & violin plots ------------------------- ###
############    -----------------------------------------    ############


### ------------------------- Main Figure X ------------------------- ###
### --------------------- Further feature plots --------------------- ###
### ------------------------ & violin plots ------------------------- ###

# ----> Output directory.
fig.markers.path <- paste0(reports.path, '/figure_on_markers')
if(!dir.exists(fig.markers.path)) dir.create(fig.markers.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4

# ---> Feature plots.
this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
tmp.data.set <- 'hlty.b'
markers.of.int <- c('FCRL5', 'FCGR2B', 'CD83', 'CD40', 'CD27', 'IGHM', 'IGHG1')
markers.of.int <- translate.ids(ids=markers.of.int, ensembl=FALSE)
tmp.file.name <- paste0(fig.markers.path, '/', obj.extended.names[tmp.data.set], '_Markers.pdf')
pdf(file=tmp.file.name)
for(marker in names(markers.of.int)){
  tmp.ggplot <- FeaturePlot(object=srt.objs.list[[tmp.data.set]], features=marker) + scale_color_gradientn(colors=this.color.scale) + labs(title=markers.of.int[marker], x='UMAP 1', y='UMAP 2')
  print(tmp.ggplot)
}
dev.off()


############    -----------------------------------------    ############
### --------------------------- Figure X ---------------------------- ###
### ------------------------ On correlations ------------------------ ###
############    -----------------------------------------    ############


### ------------------------- Main Figure X ------------------------- ###
### ------------------------ On correlations ------------------------ ###
# ----> Output directory.
fig.corr.path <- paste0(reports.path, '/figure_on_corrs')
if(!dir.exists(fig.corr.path)) dir.create(fig.corr.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4

# ---> General correlation analysis.
# @ Obtain data
# pops.per.cell.type <- list(
#   CD8=c('Teff', 'TRM', 'TCM'),
#   CD4=c('CTL', 'TRM', 'TCM')
# )
tmp.props <- lapply(X=names(srt.objs.list), FUN=function(cell.type){
  meta.data <- as.data.table(srt.objs.list[[cell.type]]@meta.data)
  # tmp.pops <- pops.per.cell.type[[cell.type]]
  tmp.lab <- clust.labs[cell.type]
  # Summary per population per donor.
  meta.data[, pop.tag:=get(tmp.lab)]
  tmp.data <- meta.data[
    !is.na(donor.id.tag),
    .(cell.count=.N),
    by=.(donor=donor.id.tag, pop.tag)
  ]
  tmp.data.2 <- meta.data[!is.na(donor.id.tag), .(total.cells=.N), by=.(donor=donor.id.tag)]
  tmp.data <- merge(x=tmp.data, y=tmp.data.2, by='donor', all=TRUE)
  tmp.data[, cell.prop:=cell.count/total.cells]
  tmp.data[, total.cells:=NULL]
  tmp.data$pop.tag <- clusts.defs[[cell.type]][tmp.data[, pop.tag]]
  return(tmp.data)
})
names(tmp.props) <- toupper(str_replace(string=names(srt.objs.list), pattern='hlty.', replacement=''))
# Further wrangling
tmp.props <- rbindlist(l=tmp.props, use.names=TRUE, idcol='cell.type')
tmp.props[, cell.count:=NULL]
tmp.props[, immune.subset:=paste(cell.type, pop.tag, sep=' ')]
tmp.props[, `:=`(cell.type=NULL, pop.tag=NULL)]
tmp.props <- spread(data=tmp.props, key=immune.subset, value=cell.prop, fill=0)

# @ Apply correlation analysis.
# corr.data.1 <- cor(all.cont.vars)
# corr.data.2 <- cor(all.cont.vars, method='spearman')
# tmp.file.name <- paste0(fig.2.path, '/CorrAnalysis.pdf')
# pdf(file=tmp.file.name)
# corrplot(corr.data.1, order="hclust") #, type = "upper")
# corrplot(corr.data.2, order="hclust") #, type = "upper")
# dev.off()

# @ Taking P-values into consideration.
corr.data.3 <- rcorr(as.matrix(tmp.props)[, 2:ncol(tmp.props)], type='spearman')
tmp.file.name <- paste0(fig.corr.path, '/CorrAnalysis_WithPVals.pdf')
pdf(file=tmp.file.name, width=10, height=10)
corrplot(corr.data.3$r, order="original", p.mat=corr.data.3$P, sig.level=0.05, insig='blank') #, type = "upper")
corrplot(corr.data.3$r, order="original", p.mat=corr.data.3$P, sig.level=0.01, insig='blank') #, type = "upper")
corrplot(corr.data.3$r, order="hclust", p.mat=corr.data.3$P, sig.level=0.05, insig='blank') #, type = "upper")
dev.off()


############    -----------------------------------------    ############
### --------------------- Figure on Correlations -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.corr.path <- paste0(reports.path, '/figure_on_corrs')
if(!dir.exists(fig.corr.path)) dir.create(fig.corr.path)


get.gene.mean.per.donor <- function(seurat.obj, gene, cluster.tag, cluster){
  # Expression data.
  exp.data <- GetAssayData(object=seurat.obj, slot='data')
  exp.data <- exp.data[gene, ]
  exp.data <- data.table(cell.bc=names(exp.data), gene=exp.data)
  # Donor metadata.
  meta.data <- as.data.table(seurat.obj@meta.data)
  meta.data[, cell.bc:=Cells(seurat.obj)]
  meta.data <- meta.data[!is.na(donor.id.tag), .(cell.bc, donor.id.tag, cluster.tag=get(cluster.tag))]
  meta.data <- meta.data[cluster.tag==cluster]
  # Merge data and calculate mean.
  tmp.data <- merge(x=meta.data, y=exp.data, by='cell.bc')
  tmp.data <- tmp.data[, .(gene=mean(gene)), by=donor.id.tag]
  colnames(tmp.data)[colnames(tmp.data)=='gene'] <- gene
  return(tmp.data)
}

tmp.data.1 <- get.gene.mean.per.donor(seurat.obj=srt.objs.list[['hlty.cd8']], gene='ENSG00000143184', cluster.tag=clust.labs['hlty.mye'], cluster='1')
tmp.data.2 <- get.gene.mean.per.donor(seurat.obj=srt.objs.list[['hlty.cd8']], gene='ENSG00000143185', cluster.tag=clust.labs['hlty.mye'], cluster='1')
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')

# tmp.ggplot <- coexpress.plot(seurat.obj=srt.objs.list[['hlty.cd8']], feature.x='ENSG00000143184', feature.y='ENSG00000143185', slot='data', filter.tags=clust.labs['hlty.cd8'], na.rm=TRUE, groups.to.filter=c('1'), keep=TRUE, feature.thold=NULL, express.thold=1, use.dp=FALSE, file.name=NULL, this.color.scale=NULL)

meta.data <- as.data.table(srt.objs.list[['hlty.mye']]@meta.data)
cdc1.enr <- meta.data[
  !is.na(donor.id.tag),
  .(
    cDC1=.SD[get(clust.labs['hlty.mye'])=='6', .N] / .N
  ),
  by=donor.id.tag
]

tmp.data <- merge(x=tmp.data, y=cdc1.enr, by='donor.id.tag')

tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=ENSG00000143184, y=ENSG00000143185)) +
  geom_point() + stat_cor() +
  labs(x='XCL1 (mean)', y='XCL2 (mean)') +
  theme_bw()

tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=ENSG00000143184, y=cDC1)) +
  geom_point() + stat_cor() +
  labs(x='XCL1 (mean)', y='% cDC1 per donor') +
  theme_bw()

tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=ENSG00000143185, y=cDC1)) +
  geom_point() + stat_cor() +
  labs(x='XCL1 (mean)', y='% cDC1 per donor') +
  theme_bw()

tmp.file.name <- paste0(fig.corr.path, '/XCL1vsXCL2_Hlty-CD8.pdf')
pdf(file=tmp.file.name)
# tmp.ggplot.2 <- tmp.ggplot + labs(x='XCL1', y='XCL2') + theme_bw()
# print(tmp.ggplot.2)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
dev.off()


############    -----------------------------------------    ############
### -------------------- Figure on summary stats -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.summ.stats.path <- paste0(reports.path, '/table_on_summ_stats')
if(!dir.exists(table.summ.stats.path)) dir.create(table.summ.stats.path)

# ---------------------------------------------------------------------->
# Name: Retrieve expression data.
# Description:
# This function will provide a dgcMatrix for a subset of cells -taken from the main environment seurat object- where the subset is taken according to a tag of interest and a specific value of such tag (i.e., the group of cells belonging to such tag value). Plus, data can be output as LogNormalized data or CPMs.
# If no more than 10 cells, the function will return NULL value.
# Arguments ------------------------->
# tag.of.interest - Character, specific tag of interest.
# group - Character, group of interest (must be a valid value for the tag of interest).
# seurat.obj and tmp.slot variable values are taken from the main environment.
# Function:
get.exp.data <- function(seurat.obj, tag.of.interest, group, data.slot='counts', norm.data='cpms'){
  # Get expression data and subset keeping only group cells.
  exp.data <- GetAssayData(object=seurat.obj, slot=data.slot, assay='RNA')
  group.cells <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.interest]==group & !is.na(seurat.obj@meta.data[, tag.of.interest])]
  if(length(group.cells)<10) return(NULL)
  exp.data <- exp.data[, group.cells]
  # Apply normalization if needed.
  if(norm.data=='cpms'){
    scale.factor <- 1000000
    size.factors <- scale.factor / Matrix::colSums(exp.data)
    exp.data <- Matrix::t(size.factors * Matrix::t(exp.data))
  }
  return(exp.data)
}
# ---------------------------------------------------------------------->

# ---> Calculate summary stats.
# Define threshold of cells to apply calculation on.
max.cells.list <- c(
  'hlty.cd4'=100000,
  'hlty.cd8'=100000,
  'hlty.nk'=30000,
  'hlty.b'=100000,
  'hlty.mye'=100000
)
# Apply calculations.
stats.per.values <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
  cat(paste0('\nPopulation: ', tmp.pop, '\n'))
  # Define inputs.
  seurat.obj <- srt.objs.list[[tmp.pop]]
  tmp.tag <- clust.labs[tmp.pop]
  # If necessary, downsample cells from the seurat objects.
  dea.max.cells <- max.cells.list[tmp.pop]
  if(length(Cells(seurat.obj)) > dea.max.cells){
    tmp.sample <- sample(x=Cells(seurat.obj), size=dea.max.cells, replace=FALSE)
    seurat.obj <- subset(x=seurat.obj, cells=tmp.sample)
  }
  # @ Tag values as groups.
  tag.values <- unique(as.character(seurat.obj@meta.data[, tmp.tag]))
  # Remove NA values.
  tag.values <- gtools::mixedsort(tag.values[!is.na(tag.values)])
  # @ Retrieve stats per group.
  stats.per.values <- lapply(X=tag.values, FUN=function(tmp.val){
    cat(paste0(tmp.val, '\n'))
    # @ Group cells' expression values.
    exp.data <- get.exp.data(seurat.obj=seurat.obj, tag.of.interest=tmp.tag, group=tmp.val)
    if(is.null(exp.data)) return(NULL)
    # @ Mean
    feat.means <- Matrix::rowMeans(exp.data)
    # @ Proportion of expressing cells.
    feat.props <- Matrix::rowSums(exp.data>0)/ncol(exp.data)
    # @ Mean of expressing cells.
    exp.data[exp.data==0] <- NA
    feat.pos.means <- Matrix::rowMeans(exp.data, na.rm=TRUE) # NA result means there's no data with value different than NA (where NA=0 in our case). Then, we need to change such values.
    feat.pos.means[is.na(feat.pos.means)] <- 0
    # @ Pack data and deliver.
    feat.stats <- data.table(feature=rownames(seurat.obj), mean=feat.means, mean.pos=feat.pos.means, prop=feat.props)
    return(feat.stats)
  })
  # Get a tag label and name accordingly.
  tag.lab <- ifelse(test=grepl(x=tmp.tag, pattern='RNA_snn_res'), yes='resolution', no=str_replace(string=tmp.tag, pattern='\\.tag', replacement=''))
  names(stats.per.values) <- paste((tag.lab), tag.values, sep='.')
  # Remove any NULL values (smallest groups with cell no < 10)
  stats.per.values <- stats.per.values[!sapply(X=stats.per.values, FUN=is.null)]
  # Condense values into a single object.
  stats.per.values <- rbindlist(l=stats.per.values, use.names=TRUE, idcol='tag.value')
  return(stats.per.values)
})
names(stats.per.values) <- names(srt.objs.list)
# Standardize format.
tmp.data <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
  # Spread data stat by stat.
  pop.data <- stats.per.values[[tmp.pop]]
  pop.data[, tag.value:=str_replace(string=tag.value, pattern='resolution.', replacement='')]
  # For mean
  tmp.data.1 <- pop.data[, .(mean=tag.value, feature, stat=mean)]
  tmp.data.1 <- spread(data=tmp.data.1, key=mean, value=stat, sep='.')
  # For mean of positive cells.
  tmp.data.2 <- pop.data[, .(pos_mean=tag.value, feature, stat=mean.pos)]
  tmp.data.2 <- spread(data=tmp.data.2, key=pos_mean, value=stat, sep='.')
  # For precentage of positive cells.
  tmp.data.3 <- pop.data[, .(prop=tag.value, feature, stat=prop)]
  tmp.data.3 <- spread(data=tmp.data.3, key=prop, value=stat, sep='.')
  # Sanity check or merge data.
  tmp.check <- (tmp.data.1[, .N] == tmp.data.2[, .N]) & (tmp.data.2[, .N] == tmp.data.3[, .N])
  if(!tmp.check) stop(paste0('Unexpected error for pop.: ', tmp.pop))
  tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='feature')
  tmp.data <- merge(x=tmp.data, y=tmp.data.3, by='feature')
  # Sanity check.
  tmp.lab <- str_replace(string=tmp.pop, pattern='tumor.|hlty.', replacement='')
  colnames(tmp.data) <- paste0(tmp.lab, '.', colnames(tmp.data))
  colnames(tmp.data)[colnames(tmp.data)==paste0(tmp.lab, '.feature')] <- 'gene'
  return(tmp.data)
})
# Merge calculations across cell types.
tmp.data <- Reduce(x=tmp.data, f=function(x, y){
  merge(x=x, y=y, by='gene', all=TRUE)
})
tmp.data.2 <- as.data.table(feature.info)
tmp.data.2 <- tmp.data.2[, .(gene=ensembl, name)]
tmp.data <- merge(x=tmp.data.2, y=tmp.data, by='gene', all.x=FALSE, all.y=TRUE)
# Report.
tmp.file.name <- paste0(table.summ.stats.path, '/SummaryStats_LungData.csv')
fwrite(file=tmp.file.name, x=tmp.data)