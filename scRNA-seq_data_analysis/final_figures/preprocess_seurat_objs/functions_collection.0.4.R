############    -----------------------------------------    ############
### --------------------- Functions collection ---------------------- ###
### -------------------- for Kusnadi et al 2021 --------------------- ###
############    -----------------------------------------    ############


# ---> About the script.
# @ Version: 0
# @ Version updates:
#   None, no stable version yet.
# @ Subversion: 1
# @ Subversion updates:
#   TBD.
# @ Updated functions:
#   None.


# 0 -------------------------------------------------------------------->
# Name: Load feature info file from cellranger output
# Description:
# Given the output folder from cellranger aggr or cellranger count (directory version), this function will load the feature info file within the directory.
# Arguments ------------------------->
# tenx.raw.data - Character, absolute path to cellranger aggr or count output directory.
# Values ---------------------------->
# Dataframe, represents the feature info table in cellranger aggr or count output.
# Function:

load.feature.info.file <- function(tenx.raw.data){
  feature.info.file <- unique(list.files(path=tenx.raw.data, pattern='features.tsv', full.names=TRUE))
  if(length(feature.info.file)!=1) stop('File with gene names -features.tsv- not properly defined.\n')
  is.zipped <- grepl(x=feature.info.file, pattern='\\.gz')
  if(is.zipped){
    # Then, read relatiosnships between feature names and ENSEMBL IDs.
    tmp.dir <- paste0(tenx.raw.data, '/tmp_with_id_', paste0(sample(x=LETTERS, size=5, replace=TRUE), collapse=''))
    tmp.cmmnd <- paste0('mkdir ',  tmp.dir, ' && cp ', feature.info.file, ' ', tmp.dir)
    system(command=tmp.cmmnd)
    tmp.file.name <- paste0(tmp.dir, '/features.tsv')
    tmp.cmmnd <- paste0('gunzip ', tmp.file.name, '.gz')
    system(tmp.cmmnd)
    feature.info <- read.delim(file=tmp.file.name, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
    tmp.cmmnd <- paste0('rm -r ', tmp.dir)
    system(tmp.cmmnd)
  }else{
    feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
  }
  return(feature.info)
}


# 0 -------------------------------------------------------------------->
# Name: Translate ID into name.
# Description:
# Given a set (vector) of ENSEMBL IDs or common gene names, this function will return the same vector with corresponding common gene names and ENSEMBL IDs, respectively, as a vector.
# ~~~* Requires the info from the data frame called feature.names (see code within seurat analysis script versions 1.4.2 or higher).
# Arguments ------------------------->
# ids - Vector with ENSEMBL IDs or common gene names as values.
# ensmebl - Logical indicatin whether ENSEMBL IDs are being input or, otherwise, common gene names. Default: TRUE
# Values ---------------------------->
# * Same vector with common gene names as vector names.
# Function:

translate.ids <- function(ids, ensembl=TRUE){
  # Define ID columns as in the dictionary.
  orig.id.col <- ifelse(test=ensembl, yes='ensembl', no='name') # Original ID.
  alt.id.col <- ifelse(test=!ensembl, yes='ensembl', no='name') # Alternative ID, to translate into.
  # Translate per ID.
  names(ids) <- sapply(X=ids, FUN=function(tmp.id){
    tmp.names <- unique(feature.info[feature.info[, orig.id.col]==tmp.id, alt.id.col])
    if(length(tmp.names)==0){
      feat.name <- NA
    }else{
      if(length(tmp.names)>1) feat.name <- paste0(tmp.names, collapse=';') else feat.name <- tmp.names
    }
    return(feat.name)
  })
  return(ids)
}


# 1 -------------------------------------------------------------------->
# Name: Dot plot.
# Of note: This function's idea was based on Seurat's output dot plot, but it was developed independently up to the previous version (see past function's file version). The size scale function implementation was copied from Seurat's function.
# Updates:
# @     Compared to the function in version 1.4, this one:
# 1. If required, applies a row-wise normalization by z-scoring.
# 2. Can set boundaries for both scales, color and size, through arguments col.min/col.max and size.min/size.max.
# 3. If required, sets a hot and cold color scale taking into account the minimum and maximum color scale values by centering at zero. This won't work if there's only negative or only positive values. (This needs further revising).
# 4. Sets a scale by taking into consideration either the size or the radius of the points in the plot. Original seurat's function does apply a radius scaling and it seems to be the best option for this kind of plots.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 4. Calculate expression mean and burst frequency for each gene across groups.
# 5. Normalize the mean expression values across the groups of interest (i.e., per gene or in a row-wise manner).
# 6. Set color and size scales by setting boundaries and plugging the appropriate arguments.
# 7. Generate a dot plot (based on seurat's style) as a ggplot, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
#     @ About data.
#         Data specificities.
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# features - Features for which to output expression values' maetrics. All should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. Valid values are the usual ones: counts, data or scale.data. Default: data.
# do.norm - Logical, indicates whether a row-wise normalization should be applied through z-scoring. Default: FALSE.
# ensembl - Logical, defines whether input features' vector (defined through argument 'fetaures') lists ENSEMBL IDs (then expecting the gene names to be defined as the vector's names). Default: FALSE.
#         Data filtering and grouping.
# groups.tag - Tag defined in the seurat object metadata whose groups will be the basis for the x axis of the product. No default.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# groups.order - Order to show the groups in the x axis. When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around. If none provided (i.e., NULL value), the function will attempt to sort them based on their names (in an alphanumeric order). Default: NULL
# filter.tag - Tag defined in the seurat object metadata whose groups will be used to filter out cells. No cell filtering is applied when this argument's set to NULL. Default: NULL.
# groups.to.filter - Character listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
#     @ About plot and output.
#         Color scale.
# this.color.scale - Character, indicates color scale to use to represent average expression. If set to NULL, a color scale that's broadly used in Vijay's lab will be used instead. If a length-1 vector that equals 'hot.and.cold', 'hot.n.cold' or 'cold.and.hot' is input, a hot and cold color scale is infered by taking into account the range of average expression values seen in the data (post normalization if required).
# col.min - Numeric, lower boundary for color scale. Must be smaller than col.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# col.max - Numeric, upper boundary for color scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Color scale.
# scale.by - Character, indicates whether a radius scaling or a standard size scaling should be applied. Accepted values are 'radius' and 'size', respectively. Default: 'radius'.
# dot.scale - Character, second value in the argument 'range' of 'scale_radius' or 'scale_size'. Default: 6.
# size.min - Numeric, lower boundary for size scale. Must be smaller than size.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# size.max - Numeric, upper boundary for size scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Output.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

dot.plot <- function(
  # @   About data.
  # Data specificities.
  seurat.obj, features, slot='data', do.norm=FALSE, ensembl=FALSE,
  # Data filtering and grouping.
  groups.tag, groups.order=NULL, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  # @   About plot and output.
  # Color scale
  this.color.scale=NULL, col.min=NULL, col.max=NULL,
  # Size scale.
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  # File.
  file.name=NULL
){
  # ---> Arguments.
  # Set color scale.
  if(is.null(this.color.scale)){
    this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  }else{
    if(this.color.scale=='hot.and.cold' | this.color.scale=='cold.and.hot' | this.color.scale=='cold.n.hot') this.color.scale <- NA
  }
  # Set size scale.
  scale.func <- switch(
    EXPR=scale.by,
    size=scale_size,
    radius=scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  # ---> Fetch data.
  # Look for all features. Make sure it's been saved to a valid slot and try to find the features there.
  if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
  feature.vals <- try(expr=as.data.table(FetchData(object=seurat.obj, vars=features, slot=slot)), silent=TRUE)
  # Kill if it was not possible to find the feature in the data slot appointed.
  if(any(class(feature.vals)=='try-error')) stop('Feature found nowhere. Please enter a valid feature.\n')
  colnames(feature.vals) <- if(ensembl) names(features) else features
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)) if(!(filter.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for filtering set to a value different than NULL, but not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)){
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag], filter.tag=seurat.obj@meta.data[, filter.tag])
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='rest']
  }
  # ---> Filter out non-desired cells.
  # According to filter tag.
  if(!is.null(filter.tag) & !is.null(groups.to.filter)){
    if(keep) tmp.data <- tmp.data[filter.tag %in% groups.to.filter] else tmp.data <- tmp.data[!filter.tag %in% groups.to.filter]
    tmp.data[, filter.tag:=NULL]
  }
  # ---> Tidy data up!
  tmp.data <- as.data.table(gather(data=tmp.data, key='feature', value='exp.val', -groups.tag))
  # ---> Apply feature thereshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[exp.val > feature.thold]
  }
  # ---> Values for color and size scales.
  tmp.data <- tmp.data[, .(exp.mean=mean(exp.val, na.rm=TRUE), burst.freq=.SD[exp.val > 0, .N]*100/.N), by=.(groups.tag, feature)]
  # ---> Row-wise normalization.
  # Apply only if required.
  if(do.norm){
    tmp.data.2 <- as.data.frame(spread(data=tmp.data[, .(groups.tag, feature, exp.mean)], key=groups.tag, value=exp.mean))
    row.names(tmp.data.2) <- tmp.data.2$feature; tmp.data.2$feature <- NULL
    tmp.data.2 <- as.data.frame(scale(x=tmp.data.2))
    tmp.data.2$feature <- row.names(tmp.data.2); tmp.data.2 <- as.data.table(tmp.data.2)
    tmp.data.2 <- gather(data=tmp.data.2, key='groups.tag', value='exp.mean', -feature)
    tmp.data <- merge(x=tmp.data[, .(groups.tag, feature, burst.freq)], y=tmp.data.2, by=c('feature', 'groups.tag'))
  }
  # ---> Set boundaries for scales and get final color definitions.
  # Color scale.
  if(!is.null(col.min) & !is.null(col.max)) if(col.min > col.max) stop('Boundaries for color scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.null(col.min)) tmp.data[exp.mean < col.min, exp.mean:=col.min]
  if(!is.null(col.max)) tmp.data[exp.mean > col.max, exp.mean:=col.max]
  if(is.na(this.color.scale)){
    scale.lims <- tmp.data[, c(max=max(exp.mean), min=min(exp.mean))]
    is.max <- names(scale.lims)[abs(scale.lims)==max(abs(scale.lims))]
    max.col.count <- 1001
    min.col.count <- floor(min(abs(scale.lims))/max(abs(scale.lims)) * max.col.count)
    max.col.scale <- heatmaply::cool_warm(max.col.count)
    min.col.scale<- heatmaply::cool_warm(min.col.count)
    if(is.max=='max'){
      this.color.scale <- c(min.col.scale[1:floor(min.col.count/2)], max.col.scale[floor(max.col.count/2):length(max.col.scale)])
    }else{
      this.color.scale <- c(max.col.scale[1:floor(max.col.count/2)], min.col.scale[floor(min.col.count/2):length(min.col.scale)])
    }
  }
  # Size scale.
  if(!is.na(size.min) & !is.na(size.max)) if(size.min > col.max) stop('Boundaries for size scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.na(size.min)) tmp.data[burst.freq<size.min, burst.freq:=size.min]
  if(!is.na(size.max)) tmp.data[burst.freq>size.max, burst.freq:=size.max]
  # ---> Set factors.
  # For groups.
  if(is.null(groups.order)){
    new.lvls <- mixedsort(x=unique(as.character(tmp.data$groups.tag)))
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
  }else{
    if(!(all(groups.order %in% unique(as.character(tmp.data$groups.tag))) & all(unique(as.character(tmp.data$groups.tag)) %in% groups.order))) stop('When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around.\n')
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=as.character(groups.order))
  }
  # For features (same order as the one provided)
  new.lvls <- if(ensembl) names(features) else features
  tmp.data$feature <- factor(x=tmp.data$feature, levels=rev(new.lvls))
  # ---> Dot plot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, size=burst.freq, fill=exp.mean)) +
    geom_point(stroke=0.4, shape=21, col='black') +
    # geom_point() +
    scale_fill_gradientn(colors=this.color.scale) +
    scale.func(range=c(0, dot.scale), limits=c(size.min, size.max)) +
    labs(x='Groups', y='Gene', fill='Average expression', size='Percent expressed') +
    theme(panel.background=element_blank(), panel.border=element_blank(), axis.line=element_line())
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}

# 2 -------------------------------------------------------------------->
# Name: Output clonotypes types values barplots across tag values.
# Description:
# Sample a group of cells from a tag value according to a established size.
# If the size is not reachable by the actual size, the function will output all of the possible cells.
# Arguments ------------------------->
# seurat.obj
# the.tag - Tag of interest (where the group belongs)
# the.group = The group of interest within the tag.
# the.size - Size for sampling.
# Function:

downsample.group <- function(seurat.obj, the.tag, the.group, the.size){
  #cat(the.group, '\n')
  potential.cells <- Cells(seurat.obj)[seurat.obj@meta.data[, the.tag]==the.group & !is.na(seurat.obj@meta.data[, the.tag])]
  # Depending if whether there are enough cells to sample.
  output.cells <- if(length(potential.cells)>the.size) sample(x=potential.cells, size=the.size) else potential.cells
  return(output.cells)
}

# 3 -------------------------------------------------------------------->
# Name: Tag analysis
# Description:
# Provide UMAP plots depicting a normalized proportion (based on each cluster) of every group of a given tag.
# Arguments ------------------------->
# seurat.obj
# tmp.tag - Tag of groups of interest.
# clusters.tag = Tag listing clusters of the seurat object.
# tag.reports.path - Absolute path where plots will be saved to.
# reports.pref - Preffix to use for output file names.
# Function:

get.tag.analysis <- function(seurat.obj, tmp.meta.data, tmp.tag, clusters.tag, tag.reports.path, reports.pref, vals.cols, output.umaps=FALSE, output.single=FALSE, downsample=TRUE, blank.complement, alpha.val=0.4){
  cat('\n\nProcess for tag', tmp.tag, '\n')

  ### ---------------------- Data preprocessing ----------------------- ###
  total.cells <- length(Cells(seurat.obj))

  # ---> Set the basic thing for any tag.
  # We get this tag's values (all except NA).
  tag.vals <- unique(as.character(seurat.obj@meta.data[, tmp.tag]))
  tag.vals <- tag.vals[!is.na(tag.vals)]
  # Amount of tag values and width length to increase based on this.
  tag.vals.no <- length(tag.vals)
  plus.width <- tag.vals.no*0.05
  # If there aren't more than one group to set comparisons, skip it.
  if(tag.vals.no < 2){ cat('Tag has a single value.'); return(NULL)}
  # Identify cells with non-NA values for this tag.
  tag.cells <- Cells(seurat.obj)[!is.na(seurat.obj@meta.data[, tmp.tag])]

  ### ------------------------ General reports ------------------------ ###
  ### ---------------------- General calculations --------------------- ###

  # ------> Tag values per cluster.
  # ---> Raw.
  tag.by.clusters <- as.matrix(table(seurat.obj@meta.data[,c(tmp.tag, clusters.tag)], exclude=c(NA, 'NA')))
  # ---> Normalized.
  # Get totals for normalization purposes.
  values.totals <- apply(X=tag.by.clusters, MARGIN=1, FUN=sum)
  values.totals <- values.totals[values.totals>0]
  tag.by.clusters <- tag.by.clusters[names(values.totals), ]
  # Identify normalization factor.
  #   Either the size of the smalles group whenever it accounts for 5% of the total or the median of the groups' sizes otherwise.
  smallest.size <- min(values.totals)
  smallest.prop <- smallest.size/sum(values.totals)
  gen.norm.factor <- ifelse(test=smallest.prop>0.05, yes=smallest.size, no=median(values.totals))
  cat('\n\n\nNormalization factor is:', gen.norm.factor, '\n\n\n')
  # Normalization (using normalization factors)
  norm.factors <- gen.norm.factor/values.totals
  tag.by.clusters.normt <- t(sapply(X=names(norm.factors), FUN=function(tag.value) tag.by.clusters[tag.value, ]*norm.factors[tag.value]))
  # Finally, add total to both tables.
  tag.by.clusters <- cbind(tag.by.clusters, values.totals)
  colnames(tag.by.clusters)[ncol(tag.by.clusters)] <- 'total'

  ### ---------------- Tag on the dimensionality plot ----------------- ###

  # ------> Define tag colors.
  if(!is.null(vals.cols)){
    tags.cols <- data.table(value=names(vals.cols), color=vals.cols)
  }

  # ------> Preprocess data.
  # Tidy table and downsampled set of cells.
  # ---> Tidy data.
  # tmp.df <- as.data.frame(tag.by.clusters.normt)
  # tmp.df$group <- rownames(tmp.df)
  # tmp.df <- gather(data=tmp.df, key='cluster', value='group.size', -group)
  # tmp.df$cluster <- factor(x=tmp.df$cluster, levels=gtools::mixedsort(unique(tmp.df$cluster)))
  # ---> Cells to depict according to nomalization-based-downsampling.
  tag.cells.set <- unlist(lapply(X=tag.vals, FUN=downsample.group, seurat.obj=seurat.obj, the.tag=tmp.tag, the.size=gen.norm.factor))

  # -----> Tag value-specific cells on dimensionality reduction plot.
  if(output.umaps) for(tmp.tag.val in tag.vals){
    cells.to.plot <- if(downsample) Cells(seurat.obj)[seurat.obj@meta.data[, tmp.tag]==tmp.tag.val & !is.na(seurat.obj@meta.data[, tmp.tag]) & Cells(seurat.obj)%in%tag.cells.set] else Cells(seurat.obj)[seurat.obj@meta.data[, tmp.tag]==tmp.tag.val & !is.na(seurat.obj@meta.data[, tmp.tag])]
    tmp.data <- cbind(seurat.obj@meta.data[Cells(seurat.obj) %in% cells.to.plot, ], seurat.obj@reductions$umap@cell.embeddings[Cells(seurat.obj) %in% cells.to.plot, ])
    # Normal version.
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=UMAP_1, y=UMAP_2)) +
      geom_point(aes_string(col=tmp.tag), size=4, alpha=alpha.val, shape=19) +
      scale_color_manual(values=vals.cols[names(vals.cols)==tmp.tag.val]) +
      geom_density_2d(color='black', size=2) +
      scale_y_continuous(expand=expansion(add=0.2)) + scale_x_continuous(expand=expansion(add=0.2)) +
      labs(title=tmp.tag.val, x='UMAP 1', y='UMAP 2', col='Cluster')
    # Output.
    tmp.pref <- paste0(reports.pref, '_TagVal-', tmp.tag.val)
    publish.plot(tmp.ggplot=tmp.ggplot, output.path=tag.reports.path, file.name=tmp.pref, type='tiff', blank.comp=blank.complement, do.legend=FALSE, height=5, width=5)
  }
  if(output.single){
    cells.to.plot <- Cells(seurat.obj)[Cells(seurat.obj)%in%tag.cells.set]
    tmp.data <- cbind(seurat.obj@meta.data[Cells(seurat.obj) %in% cells.to.plot, ], seurat.obj@reductions$umap@cell.embeddings[Cells(seurat.obj) %in% cells.to.plot, ])
    # Normal version.
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=UMAP_1, y=UMAP_2)) +
      geom_point(aes_string(col=tmp.tag), size=4, alpha=alpha.val, shape=19) +
      scale_color_manual(values=vals.cols) +
      geom_density_2d(color='black', size=2) +
      scale_y_continuous(expand=expansion(add=0.2)) + scale_x_continuous(expand=expansion(add=0.2)) +
      labs(x='UMAP 1', y='UMAP 2', col='Group')
    # Output.
    tmp.pref <- paste0(reports.pref, '_AllTagVals')
    publish.plot(tmp.ggplot=tmp.ggplot, output.path=tag.reports.path, file.name=tmp.pref, type='tiff', blank.comp=blank.complement, do.legend=FALSE, height=5, width=5)
  }
  # -----> Return table of normalized values.
  return(tag.by.clusters.normt)
}

# 4 -------------------------------------------------------------------->
# Name: Heatmap with DEA results sorted according to sets' values and sizes.
# Description:
# Given a table (provided either as a data frame or a data table) listing the metrics ussually gotten when applying multiple pairwise comparisons through DEA (mainly, adjusted p-value, LFC and mean values for the groups that were compared), this function will provide a tidied version of such table listing the DEGs in an ordered way that depends on both the reference elements of each comparison (set value) and the size of groups of elements that share a DEG (set size). When a heatmap is plotted based on this order (provided an order of the elements), it must have a scalonated pattern. In order to get such results, the function will:
# 1. Tidy the data up keeping relevan columns and calculate the size set each DEG belongs to.
# 2. Set order, first according to set size, then to LFC and finally to elements.
# Arguments ------------------------->
# dea.results - DEA results table. Data frame or data table, must contain all DEA metrics as well as some identifier for each DEG. No default.
# feat.id - Character, name of the column listing DEG's ID. Default: 'gene.id'
# element - Character, name of the column listing each element value (i.e., the value for which the DEA was applied (either as reference or not)). Make sure this column's name is different than 'element', since the function will not work otherwise. No default.
# vals.order - Character vector containing the order for the element values desired to show on the heatmap. No default.
# p.adj - Character, name of the column listing DEG's adjusted p values. Default: 'p.adj'
# lfc - Character, name of the column listing DEG's LFC. Default: 'lfc'
# mean.tag - Character, pattern to identify columns listing DEA group means. No default.
# element.eq.means - Logical, indicates whether element values are the same values for which a mean is provided. If TRUE, the function will try to get rid of the information taht's repeated along the process. Default: FALSE.
# file.name - To implement. If not NULL, the function should output the heatmap, too.
# Value ----------------------------->
# Datatable listing DEGs main metrics ordered according to both element values and set sizes.
# Function:

dea.heatmap <- function(dea.results, feat.id='gene.id', element, vals.order, p.adj='p.adj', lfc='lfc', mean.tag, element.eq.means=FALSE, file.name=NULL){
  # ---> Preflights.
  if(!is.data.table(dea.results)) dea.results <- as.data.table(dea.results)
  # ---> Define stuff
  # Stats.
  mean.cols <- colnames(dea.results)[str_detect(string=colnames(dea.results), pattern=mean.tag)]
  stat.cols <- c(p.adj, lfc, mean.cols)
  col.vals <- str_replace(string=mean.cols, pattern=mean.tag, replacemen='')
  # ---> Tidy data up!
  stats.vals <- apply(X=as.matrix(dea.results[, ..stat.cols]), MARGIN=1, FUN=paste, collapse=';')
  dea.results[, stats:=stats.vals]
  tmp.data <- dea.results[, .(gene.id=get(feat.id), element=get(element), stats)]
  tmp.data <- spread(data=tmp.data, key='element', value='stats', fill=paste('0;0', paste0(rep(x='0', times=length(mean.cols)), collapse=';'), sep=';'))
  for(element.val in vals.order){
    tmp.data <- separate(data=tmp.data, col=element.val, into=paste(c('p.val', 'lfc', paste0('mean.', col.vals)), element.val, sep='.'), sep=';', convert=TRUE)
  }
  # Discard repeated information when the mean groups (values) are the same as the element values.
  if(element.eq.means){
    mean.cols <- paste('mean', col.vals, col.vals, sep='.')
    if(!all(mean.cols %in% colnames(tmp.data))) stop('It was requested to remove repeated information, though mean column values don\'t seem to be the same as the element values (could be that not all element values have a mean column). Please set element.eq.means to FALSE and see if your results are still favorable. Otherwise, check your input.\n')
    other.cols <- colnames(tmp.data)[!str_detect(string=colnames(tmp.data), pattern='mean.')]
    cols.to.keep <- c(other.cols, mean.cols)
    tmp.data <- tmp.data[, ..cols.to.keep]
    mean.cols <- str_replace(string=mean.cols, pattern=paste0('.', col.vals, '$'), replacement='')
    cols.to.rep <- c(other.cols, mean.cols)
    colnames(tmp.data) <- cols.to.rep
  }
  # Identify for what cell types a gene is a DEG as well as the total amount of cell types it is DEG for.
  lfc.cols <- grep(x=colnames(tmp.data), pattern='lfc', value=TRUE)
  p.cols <- grep(x=colnames(tmp.data), pattern='p.val', value=TRUE)
  rows.to.keep.1 <- apply(X=tmp.data[, ..lfc.cols], MARGIN=1, FUN=function(row) return(row>0.25))
  rows.to.keep.2 <- apply(X=tmp.data[, ..p.cols], MARGIN=1, FUN=function(row) return(row<0.05))
  rows.to.keep.mat <- rows.to.keep.1&rows.to.keep.2
  row.names(rows.to.keep.mat) <- str_replace(string=row.names(rows.to.keep.mat), pattern='lfc\\.', replacement='')
  # While doing so, filter for genes that are significantly differentially expressed for at least one cell type (apply to both objects).
  rows.to.keep <- apply(X=rows.to.keep.mat, MARGIN=2, FUN=any)
  tmp.data <- tmp.data[rows.to.keep]
  rows.to.keep.mat <- rows.to.keep.mat[, rows.to.keep]
  # Add sets sizes and values.
  degs.sets.info <- t(apply(X=rows.to.keep.mat, MARGIN=2, FUN=function(col){
    tmp.vals <- row.names(rows.to.keep.mat)[col]
    tmp.size <- length(tmp.vals)
    tmp.vals <- paste0(sort(tmp.vals), collapse=';')
    return(c(tmp.vals, tmp.size))
  }))
  colnames(degs.sets.info) <- c('set.values', 'set.size')
  tmp.data <- cbind(tmp.data, degs.sets.info)
  # ---> Set order.
  # Get sets' sizes and values to apply iterations.
  set.sizes <- 1:(dea.results[, length(unique(get(element)))])
  set.vals <- dea.results[, unique(get(element))]
  if(!all(sort(set.vals)==sort(vals.order))) stop('Not all element values were provided in argument \'vars.order\'')
  set.vals <- vals.order
  vals.order <- 1:length(set.vals)
  names(vals.order) <- set.vals
  # Get DEGs' order.
  set.degs <- lapply(X=set.sizes, FUN=function(tmp.size){
    # Get DEGs set for the size of a given iteration.
    set.degs <- tmp.data[set.size==tmp.size]
    # Then, identify possible values for the whole size set.
    sets.vals <- str_split(string=set.degs[, set.values], pattern=';', simplify=TRUE)
    # Order size set according to values.
    set.degs <- lapply(X=set.vals, FUN=function(tmp.val){
      # Idenitfy values that have already been taken into account in previos iterations.
      to.exclude <- set.sizes[set.vals == tmp.val]-1
      if(to.exclude>0) to.exclude <- 1:to.exclude
      to.exclude <- set.vals[to.exclude]
      # Get DEGs size subset for the value of a given iteration, disregarding the DEGs belonging to any value to exclude.
      rows.to.keep <- apply(X=sets.vals, MARGIN=1, FUN=function(row) return(tmp.val %in% row & !any(to.exclude %in% row)))
      set.degs <- set.degs[rows.to.keep]
      sets.vals <- sets.vals[rows.to.keep, ]
      # Check we've got info to work with.
      if(all(!rows.to.keep)) return(NA)
      # Identify LFC columns to consider to get a mean LFC value.
      iter.set.vals <- set.degs[, unique(set.values)]
      for(tmp.val.set in iter.set.vals){
        tmp.lfc.cols <- paste0('lfc.', str_split(string=tmp.val.set, pattern=';', simplify=TRUE))
        set.degs[, tmp.avg.lfc:=rowMeans(set.degs[, ..tmp.lfc.cols])]
        set.degs[set.values==tmp.val.set, avg.lfc:=tmp.avg.lfc]
      }
      set.degs[, tmp.avg.lfc:=NULL]
      # Check we've got dimensions to work over.
      if(is.null(dim(sets.vals))) return(set.degs)
      # Set order according to mean LFC while keeping the order for object having set's values, too.
      set.degs[, idx:=1:.N]
      setorderv(x=set.degs, cols='avg.lfc', order=-1)
      sets.vals <- sets.vals[set.degs[, idx], ]
      set.degs[, idx:=NULL]
      # Set order according to set values.
      sets.vals <- as.data.table(t(apply(X=sets.vals, MARGIN=1, FUN=function(row) return(sort(vals.order[row])))))
      row.names(sets.vals) <- as.character(1:nrow(sets.vals))
      colnames(sets.vals) <- as.character(1:ncol(sets.vals)); cols.to.order <- as.character(1:ncol(sets.vals)); sets.vals[, idx:=1:.N]
      for(tmp.col in rev(cols.to.order)) setorderv(x=sets.vals, col=tmp.col, order=1)
      set.degs <- set.degs[sets.vals[, idx]]
      # return
      return(set.degs)
    })
    set.degs <- set.degs[!is.na(set.degs)]
    set.degs <- rbindlist(l=set.degs, use.names=TRUE, idcol=FALSE)
    return(set.degs)
  })
  set.degs <- rbindlist(l=set.degs, use.names=TRUE, idcol=FALSE)
  # ---> Signature depicted by a heatmap.
  # To implement.
  # Return sorted data table.
  return(set.degs)
}


# 5 -------------------------------------------------------------------->
# Name: Get density values.
# Taken from Ciro's clever functions.
# Description:
# Given the x values and the y values (compraising a grid containing 'n' cells -see parameter 'n'-) for a set of 2D points, this function will calculate the density values of such a grid and then will define the 2D density values for each individual point.
# Arguments ------------------------->
# x - x values for the set of points.
# y - y values for the set of points.
# ... - Other arguments taken by function 'kde2d' from the package 'MASS'. Among them, check parameter 'n' which basically defines the number of cells to be contained in the grid (see https://rdrr.io/cran/MASS/man/kde2d.html).
# Value ----------------------------->
# Matrix with 2D density values (columns) for each input point (rows).
# Function:

get.densities <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# 6 -------------------------------------------------------------------->
# Name: Co-expression plot.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 3. Define expression-related clasess for each cell by checking whether expression for each feature is larger than a set threshold (express.thold). Possible classes for each cell are: Double positive (DP, expresses both features over the threshold); Double negative (DN, expresses none of the features over the threshold); and single positive for either feature (SP.X or SP.Y when expressing the feature depicted in the x coordinate or the y coordinate of the scatter plot, respectively).
# 4. Calculate the fraction of each expression-related class to depict on the final plot.
# 5. Calculate the 2D kernel density for each cell.
# 6. Generate a ggplot object depicting a scatter plot with density as color, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# feature.x & feature.y - Features for which to depict their expression values in scatter plot (x and y axes, respectively). Both should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. NULL indicates 'data', 'cpm' indicates log2(CPM + 1) relative counts, while the rest are the usual ones: counts, data or scale.data. No other value is allowed. Default: data.
# filter.tags - Tags defined in the seurat onject metadata whose groups will be used to filter out cells. No cell filtering is applied when it's set to NULL. Default: NULL.
# na.rm - Logical, indicates whether NA values defined the filter.tag column should be discarded for the tag listing groups. Default: TRUE.
# groups.to.filter - Character or list, listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. When more than one tag are provided for filtering, this argument must be a list with the same dimensions as that argument. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
# express.thold - Threshold used to define expression-related clasess as described above. Default: 1
# use.dp = Logical, indicates whether only the x and y values of cells classified as double positives (DP) should be used to calculate 2d kernel density values. If TRUE, the density value for the cells belonging to any other class is set to 0. Default: FALSE.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# this.color.scale - Character. Color scale or color values to use. To be implemented, i.e., is totally ignored for now.
# Dependencies ---------------------->
# Function: get.densities to calculate 2D kernel densities for each cell.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# To do ----------------------------->
# Consider a specific color scale.
# Function:

coexpress.plot <- function(seurat.obj, feature.x, feature.y, slot='data', filter.tags=NULL, na.rm=TRUE, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, express.thold=1, use.dp=FALSE, file.name=NULL, this.color.scale=NULL){
  # ---> Arguments and checks.
  if(is.null(this.color.scale)) this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # To implement.
  if(is.null(slot)) slot <- 'data'
  both.feats <- c(feature.x, feature.y)
  # Define CPMs-related arguments..
  if(slot == 'cpm'){
    slot <- 'counts'
    cpm <- TRUE
  }else{
    cpm <- FALSE
  }
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # Feature thresholds (for cells filtering). Get a 2-length vector if a single, non-NULL value provided.
  if(!is.null(feature.thold)){
    if(length(feature.thold)==0 & length(feature.thold)>2) stop('Feature threshold must have length > 0 and < 3.\n')
    if(length(feature.thold)==1) feature.thold <- c(feature.thold, feature.thold)
  }
  # Expression-related threshold. Must be a single value (positive scalar).
  if(length(express.thold)!=1 | express.thold<=0) stop('Expression threshold (argument \'express.thold\') value must be a single value (positive scalar).\n')
  # ---> Fetch data.
  # Look for feature.
  # Make sure we've got a valid slot and try to find the feature there.
  if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
  # Then, make sure we've got both features defined in whatever data slot that was defined.
  tmp.assess <- both.feats %in% rownames(GetAssayData(object=seurat.obj, slot=slot, assay='RNA'))
  if(!all(tmp.assess)){
    tmp.assess <- both.feats[!tmp.assess]
    stop(paste0('Not all features defined in slot requested, ', slot, ', as follows:\n', paste0(tmp.assess, collapse='\n'), '\n'))
  }
  # Proceed if all good.
  if(cpm){
    feature.vals <- RelativeCounts(data=GetAssayData(object=seurat.obj, slot=slot, assay='RNA'), scale=1e6, verbose=FALSE)
    feature.vals <- log2(t(as.matrix(feature.vals[both.feats, ])) + 1)
    feature.vals <- as.data.table(feature.vals)
  }else{
    feature.vals <- as.data.table(FetchData(object=seurat.obj, vars=both.feats, slot=slot))
  }
  colnames(feature.vals) <- c('feature.x', 'feature.y')
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
  }
  # If necessary, join altogether.
  if(!is.null(filter.tags)) tmp.data <- cbind(tmp.data, feature.vals) else tmp.data <- feature.vals
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Apply feature threshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[feature.x > feature.thold[1]]
    tmp.data <- tmp.data[feature.y > feature.thold[2]]
  }
  # ---> Define expression-related classes.
  # @  Define classes in general data,
  # Double positives.
  tmp.data[feature.x>express.thold & feature.y>express.thold, express.class:='DP']
  # Single positives (either for x or y).
  tmp.data[feature.x>express.thold & feature.y<=express.thold, express.class:='SP.X']
  tmp.data[feature.x<=express.thold & feature.y>express.thold, express.class:='SP.Y']
  # Double negatives.
  tmp.data[feature.x<=express.thold & feature.y<=express.thold, express.class:='DN']
  # @  Calculate fraction of each class.
  class.counts <- tmp.data[, .(class.fraction=.N/tmp.data[, .N]), by=express.class]
  # @  Get caption on classess' fractions.
  tmp.caption <- class.counts[, paste(express.class, round(x=class.fraction, digits=2), sep=' - ', collapse='; ')]
  # ---> Calculate densities.
  # Depends on whether all cells should be used or only double positive cells.
  if(use.dp){
    dens.no <- min(c(tmp.data[express.class=='DP', .N], 100))
    tmp.data[express.class=='DP', density:=get.densities(x=tmp.data[express.class=='DP', feature.x], y=tmp.data[express.class=='DP', feature.y], n=dens.no)]
    tmp.data[is.na(density), density:=0]
  }else{
    dens.no <- min(c(tmp.data[, .N], 100))
    tmp.data[, density:=get.densities(x=tmp.data[, feature.x], y=tmp.data[, feature.y], n=dens.no)]
  }
  # ---> Scatter plot depicting density, a.k.a. Co-expression plot.
  # Get ggplot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=feature.x, y=feature.y)) + geom_density2d(colour = "#c0c5ce") + geom_point(aes(col=density), size=0.2) + viridis::scale_color_viridis(option = "magma") + labs(x=feature.x, y=feature.y, col='Density', caption=tmp.caption)
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}

# 7 -------------------------------------------------------------------->
# Name: Module scoring.
# Dependencies:
#   Seurat and stringr, ggplot2, stringr, data.table, R_handy_functions
# Functions dependencies:
#   translate.ids, among others.
# Description:
# You can find an extended explanation on the AddModuleScore function from Seurat in the script associated to this module function.
# This program takes a seurat object and a directory where a set of files listing gene signatures should be and calculates module score for each signature.
# Arguments ------------------------->
# seurat.obj - seurat object.
# module.feats.dir - String for the absolute path to the directory storing one or more gene-signature lists. Each file should have a csv format (one column which must be named 'feature') with rows listing the set of features (usually gene names) that define the module and should be part of the seurat object; if they're not, they will be filtered out. Each file name must follow the next pattern '<Signature name -spaced substituted by underscores->_signature.csv
# signal.thold <- Signal threshold to pick positive and negative cells for the signal (i.e., to define module tag). Value ignored if flag FPR is activated (i.e., if it's a value different than NULL).
# is.ensembl - Logical indicating whether input seurat object has ENSEMBL IDs defined as gene IDs instead of common genen names.
# reports.path - Absolute path to where filtered gene signatures should be stored.
# Value ----------------------------->
# Seurat object with module scores calculated.
# Function:

get.module.scores <- function(seurat.obj, module.feats.dir, signal.thold, is.ensembl=FALSE, reports.path=NULL){
  cat('\n\n')
  cat('############    -----------   Module Scoring    -----------    ############\n')
  cat('############    -Based on module-defining genes expression-    ############\n')

  ### --------------------------- Arguments --------------------------- ###
  this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.

  cat('\n\n')
  ### --------------------------- Load data --------------------------- ###
  cat('### --------------------------- Load data --------------------------- ###\n')
  # Module-defining genes file.
  module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
  modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
    tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
    if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
    tmp.feats <- tmp.feats$feature
    return(tmp.feats)
  })
  names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
  names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
  names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
  modules.names <- names(modules.feats)
  cat('Module-defining features loaded...\n')

  ### ----------------------- Module definition ----------------------- ###
  cat('### ----------------------- Module definition ----------------------- ###\n')
  data.path <- if(!is.null(reports.path)) paste0(reports.path, '/data') else NULL
  if(!is.null(data.path)) if(!dir.exists(data.path)) dir.create(data.path)

  # ---> Module-defining genes.
  # Make sure that all module-defining genes are part of the seurat object.
  cat('Filtering non-defined features...\n')
  modules.feats <- lapply(X=names(modules.feats), FUN=function(tmp.module){
    module.feats <- modules.feats[[tmp.module]]
    # Translate IDs if necessary.
    if(is.ensembl) module.feats <- translate.ids(ids=module.feats, ensembl=FALSE) else names(module.feats) <- module.feats
    # Define IDs as a data frame.
    module.feats <- module.feats[!is.na(names(module.feats))]
    module.feats <- data.frame(features=names(module.feats), original.features=module.feats, present=(names(module.feats) %in% rownames(seurat.obj)), stringsAsFactors=FALSE)
    module.feats <- module.feats[!is.na(module.feats$features), ]
    # Make sure there are any features present to continue.
    tmp.check <- all(!module.feats$present)
    if(tmp.check){
      tmp.warning <- paste0('There are not enough features from module ', tmp.module, '. It will be skipped.\n')
      warning(tmp.warning)
      return(NA)
    }
    # Then, output final list of module-defining genes.
    if(!is.null(data.path)){
      tmp.file.name <- paste0(data.path, '/', str_replace_all(string=toupper(tmp.module), pattern='\\.', replacement='_'), '_ModuleDefiningGenesInExpressionData.csv')
      write.csv(x=module.feats, file=tmp.file.name, quote=FALSE, row.names=FALSE)
    }
    # Check at least 90% of the original siganture features are present in the seurat object.
    if(sum(module.feats$present)/nrow(module.feats) < 0.85) warning(paste0('\nFor module/signature ', tmp.module, ', we couldn\'t recover for the seurat object more than 85% of genes in the original list. Something may be wrong either with the list or the seurat object. Please be careful about that when interpreting the results.\n'))
    # Output final features.
    return(module.feats$features[module.feats$present])
  })
  modules.names <- modules.names[!is.na(modules.feats)]
  modules.feats <- modules.feats[!is.na(modules.feats)]
  names(modules.feats) <- modules.names
  cat('List of module-defining features depicting presence in the gene expression data output.\n')

  cat('\n\n')
  ### -------------------- Module score inference --------------------- ###
  cat('### -------------------- Module score inference --------------------- ###\n')
  # ---> Module score.
  # Add module score across cells.
  seurat.obj <- AddModuleScore(object=seurat.obj, features=modules.feats, name=modules.names)
  # Fix module names given by seurat by removing the suffix number.
  modules.names <- paste0(modules.names, '.score')
  names(modules.names) <- paste0(str_replace(string=modules.names, pattern='\\.score', replacement=''), 1:length(modules.names))
  for(tmp.module in names(modules.names)){
    colnames(seurat.obj@meta.data) <- str_replace(string=colnames(seurat.obj@meta.data), pattern=tmp.module, replacement=modules.names[tmp.module])
  }
  names(modules.names) <- str_replace(string=names(modules.names), pattern='\\d+$', replacement='')
  names(modules.names) <- str_replace_all(string=names(modules.names), pattern='\\.', replacement='_')

  return(seurat.obj)
}

# 8 -------------------------------------------------------------------->
# Taken from the in-dev version of the functions script.
# Name: Violin plot.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tags set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of one or multiple tags (filter.tags). See details below as to how to indicate filtring based on multiple tags.
# 4. Set a color scale for each group (either the final identities -"color"- or a stat from mean, median or percentage of expressing cells).
# 5. Generate a ggplot, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# feature - Feature for which to output density. Should be defined either in the metadata (slot could either be NULL for straight and unambiguous picking or simply not be defined in any data slot) or as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. NULL indicates 'data', 'cpm' indicates log2(CPM + 1) relative counts, while the rest are the usual ones: counts, data or scale.data. No other value is allowed. Default: data.
# groups.tag - Tag defined in the seurat onject metadata whose groups will be the basis for the x axis of the product. No default.
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# groups.label - Character. When a non-NULL values provided, indicates that the groups of interest provided throug argument 'groups.of.int' should be in turn merged into a single group and the value given through this argument will be used as a label for such a single group. Otherwise, a NULL value indicates that the groups of interest shouldn't be merged. Default: NULL.
# filter.tags - Chracter vector, which should list the tag or tags in the seurat object metadata whose groups will be used to filter out cells. When filtering should be applied based on multiple tags, their groups to be considered to keep or exclude cells must be passed as chracter vectors embedded within a list; then, this character vector must have the same length as that of the list provided through the argument 'groups.to.filter'. No cell filtering is applied when it's set to NULL. Default: NULL.
# groups.to.filter - List of vectors for the set of tag groups that'll be kept ot removed according to the argument 'keep' value. Must have the same length as the vector provided in 'filter.tags'. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
# color - Value for color scale of violin, either the final identities -"color"- or a stat from mean -"mean"-, median -"median"-, burst frequency (percentage of expressing cells) -"burst.frequency"- or burst size (mean of expressing cells, where a cell is considered as such when expression > threshold, (see argument size.thold)) -"burst.size"-. Default: 'color'.
# size.thold - Threshold over which cells are going to be considered as 'expressing or not'. Ignored when color different than 'burst.size'. Default: 0
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# adjust.val - Value for argument 'adjust' of geom_violin. Default: 0.8. See details at https://ggplot2.tidyverse.org/reference/geom_violin.html
# trim.val - Value for argument 'trim' of geom_violin. Default: FALSE. See details at https://ggplot2.tidyverse.org/reference/geom_violin.html
# plot.limits - Value for argument adjust of scale_fill_gradientn. Only valid if argument color's value is different than 'color'. Default: NULL. See details at: https://ggplot2.tidyverse.org/reference/scale_gradient.html
# add.points - Logical, indicates whether jittered points should be depicted or not in the violin plot. If TRUE, only 20% of the data points (picked randomly) will be depicted. Default: FALSE.
# this.color.scale - Character. Color scale or color values to use.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

vln.plot <- function(seurat.obj, feature, slot='data', groups.tag, na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, vln.type='violin', color='color', size.thold=0, file.name=NULL, adjust.val=0.8, trim.val=FALSE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL, do.fill=TRUE){
  # ---> Arguments and checks.
  # Argument defining violin plots' fill color.
  if(is.null(this.color.scale) & color!='color') this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.
  # Define CPMs-related arguments.
  if(is.null(slot)){
    cpm <- FALSE
  }else{
    if(slot == 'cpm'){
      slot <- 'counts'
      cpm <- TRUE
    }else{
      cpm <- FALSE
    }
  }
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # ---> Fetch data.
  # Look for feature.
  # First, look for it in the metadata if slot is set to NULL.
  if(is.null(slot)){
    if(!(feature %in% colnames(seurat.obj@meta.data))) stop('Slot set to NULL, but feature required was not found in seurat object\'s metadata. Were you trying to look for it any slot?\n')
    feature.vals <- data.table(feature=seurat.obj@meta.data[, c(feature)])
  }else{
    # If set to a value different than NULL, make sure it's a valid slot and try to find the feature there.
    if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
    # Proceed if all good, though process branches based on CPM possible requirement.
    if(cpm){
      if(!feature %in% row.names(seurat.obj)) stop('Feature could not be found in the count slor, which is required to calculate CPMs as requested. May want to try a different slot?\n')
      feature.vals <- RelativeCounts(data=GetAssayData(object=seurat.obj, slot=slot, assay='RNA'), scale=1e6, verbose=FALSE)
      feature.vals <- log2(feature.vals[feature, ] + 1)
      feature.vals <- as.data.table(feature.vals)
      colnames(feature.vals) <- 'feature'
    }else{
      feature.vals <- try(expr=as.data.table(FetchData(object=seurat.obj, vars=feature, slot=slot)), silent=TRUE)
      # Only if it was not possible to find the feature in the data slot appointed, try again in the metadata.
      if(any(class(feature.vals)=='try-error')){
        if(!(feature %in% colnames(seurat.obj@meta.data))) stop('Feature found nowhere. Please enter a valid feature. Alternatively, package \'data.table\' may not be imported.\n')
        feature.vals <- data.table(feature=seurat.obj@meta.data[, c(feature)])
      }else{
        # If not, set the appropriate column name for format consistency.
        colnames(feature.vals) <- 'feature'
      }
    }
  }
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
    if(color %in% colnames(seurat.obj@meta.data)){
      tmp.data <- cbind(data.table(groups.tag=seurat.obj@meta.data[, groups.tag], color=seurat.obj@meta.data[, color]), tmp.data)
    }else{
      tmp.data <- cbind(data.table(groups.tag=seurat.obj@meta.data[, groups.tag]), tmp.data)
    }
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='Rest']
  # Similarly, merge groups of interest if their label is provided.
    if(!is.null(groups.label)){
      tmp.data[groups.tag %in% groups.of.int, groups.tag:=as.character(groups.label)]
      # Set factors leaving 'Rest' label right at the end of the levels' order.
      new.lvls <- c(as.character(groups.label), 'Rest')
      tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
    }
  }
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Apply feature thereshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[feature > feature.thold]
  }
  # ---> Color scale.
  if(color=='color'){ # Color per group.
    tmp.data[, color:=groups.tag]
  }else{
    if(!'color' %in% colnames(tmp.data)){ # Stats
      switch(EXPR=color,
        "median"={
          tmp.dt <- tmp.data[, .(color=median(feature, na.rm=TRUE)), by=groups.tag]
        },
        "mean"={
          tmp.dt <- tmp.data[, .(color=mean(feature, na.rm=TRUE)), by=groups.tag]
        },
        "burst.frequency"={ # Percent of expresing cells.
          tmp.dt <- tmp.data[, .(color=.SD[feature > 0, .N]*100/.N), by=groups.tag]
        },
        "burst.size"={ # Mean of expressing cells.
          tmp.dt <- tmp.data[feature>size.thold, .(color=mean(feature, na.rm=TRUE)), by=groups.tag]
        },
        stop(paste0('Valid value for "color" must be provided among:\n\t', paste0(c('median', 'mean', 'burst.frequency', 'burst.size', 'color'), collapse=', '), '\n'))
      )
      tmp.data <- merge(x=tmp.data, y=tmp.dt, by='groups.tag')
    }
  }
  # ---> Points
  # 20% of the points per group are going to be added.
  if(add.points){
    rows.to.keep <- sample(x=1:nrow(tmp.data), size=nrow(tmp.data)*0.05)
    tmp.data.sample <- tmp.data[rows.to.keep, ]
  }
  # ---> Violin plot.
  # Get ggplot.
  switch(
    EXPR=vln.type,
    'violin'={
      tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, fill=color)) + geom_violin(alpha=0.8, trim=trim.val, adjust=adjust.val) + geom_boxplot(width=0.05, alpha=0.7, outlier.shape=NA) + labs(x=paste0('Groups (', groups.tag, ')'), y=feature)
      # Adjust color scale.
      if(!is.null(this.color.scale) & color=='color') tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=this.color.scale)
      if(color!='color') tmp.ggplot <- tmp.ggplot + scale_fill_gradientn(colors=this.color.scale, limits=plot.limits)
      # Add points if necessary.
      if(add.points) tmp.ggplot <- tmp.ggplot + geom_jitter(data=tmp.data.sample, color='black', size=0.2, alpha=0.6, width=0.2)
    },
    'density'={
      if(do.fill){
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=feature, fill=groups.tag)) + geom_density(alpha=0.8) +
        labs(x=feature, y='Density', fill=paste0('Groups (', groups.tag, ')'))
      }else{
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=feature, col=groups.tag)) + geom_density(alpha=0.8, size=1.5) +
        labs(x=feature, y='Density', color=paste0('Groups (', groups.tag, ')'))
      }
      # Adjust color scale.
      if(!is.null(this.color.scale)){
        # Confirm that entries make sense. Else, output a warning and proceed with default ggplot colors.
        tmp.check <- tmp.data[, all(groups.tag %in% names(this.color.scale))]
        if(!tmp.check){
          warning('The names of the attempted values to color density plot do not match 1:1 with the values of the tag of interest of the metadata.\n')
        }else{
          if(do.fill){
            tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=this.color.scale)
          }else{
            tmp.ggplot <- tmp.ggplot + scale_color_manual(values=this.color.scale)
          }
        }
      }
    },
    'split'={
      tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, fill=color)) + geom_split_violin() + labs(x=paste0('Groups (', groups.tag, ')'), y=feature, fill='')
      # Adjust color scale.
      # if(!is.null(this.color.scale) & color=='color') tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=this.color.scale)
      # if(color!='color') tmp.ggplot <- tmp.ggplot + scale_fill_gradientn(colors=this.color.scale, limits=plot.limits)
      # # Add points if necessary.
      # if(add.points) tmp.ggplot <- tmp.ggplot + geom_jitter(data=tmp.data.sample, color='black', size=0.2, alpha=0.6, width=0.2)
    },
    stop('No appropriate violin-like geom was defined. Valid entires are: "violin," "density" and "split".\n')
  )
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}

# 9 -------------------------------------------------------------------->
# Name: Include PS information.
# Description:
# Given a data table describing the metadata of a dataset that includes clonotype, donor class and antigen reactivity information under the columns 'clonotype.tag', 'donor.id.tag' and 'virus.tag' , respectively, this function will include the PS information for each cell in such an object.
# Arguments ------------------------->
# meta.data - Data table, with the specifications mentioned in this function's description.
# Function:

include.ps.info <- function(meta.data){
  # Get PS info.
  ps.info <- meta.data[
    !is.na(clonotype.tag) & !is.na(donor.id.tag),
    .(
      ag.size=.N
    ),
    by=.(
      ag=virus.tag,
      clone.id=clonotype.tag,
      donor.id=donor.id.tag
    )
  ]
  ps.info <- ps.info[
    ,
    .(
      pr.tag=ifelse(test=sum(ag.size>=2)>1, yes='pR', no='non-pR')
    ),
    by=.(
      clone.id,
      donor.id
    )
  ]
  # Include in current metadata.
  meta.data <- merge(x=meta.data, y=ps.info, by.x=c('clonotype.tag', 'donor.id.tag'), by.y=c('clone.id', 'donor.id'), all.x=TRUE, all.y=FALSE, sort=FALSE)
  return(meta.data)
}


# 10 ------------------------------------------------------------------->
# Name: Replace HTO 'NA' values with the classes provided by Demuxlet
# Description:
# Replace HTO 'NA' values (for those cells without HTO information) with the classes provided by Demuxlet.
# Arguments ------------------------->
# seurat.obj - Seurat object.
# set.lab
# tissue.type
# aggr.pffx.lab
# Function:

rep.hto.negatives.with.demuxlet <- function(seurat.obj, set.lab, tissue.type='Normal-lung', aggr.pffx.lab){
  # Predefinitions.
  cell.type.translate <- c( # Labels given to aggr table files.
    cd4='CD4',
    cd8='CD8',
    mye='Myeloid',
    nk='NKAndB',
    b='NKAndB'
  )
  # @ Demuxlet-related suffixes
  demuxlet.sffxs <- c(
    CD4='01-11-2023_alpha-0.5_DS-10e6',
    CD8='01-11-2023_alpha-0.5_DS-10e6',
    Myeloid='01-11-2023_alpha-0.5_DS-10e6',
    NKAndB='01-11-2023_alpha-0.5_DS-10e6'
  )
  # Determine cell type and complete aggr table ID.
  cell.type <- cell.type.translate[str_extract(string=set.lab, pattern='cd4$|cd8$|mye$|b$|nk$')]
  aggr.id <- paste0(aggr.pffx.lab, cell.type, '_', tissue.type)
  demuxlet.sffx <- demuxlet.sffxs[cell.type]
  # Read aggr table and run.
  aggr.table <- fread(file=aggr.files[aggr.id])
  # Modify aggr table according to identify files output by Demuxlet.
  seq.pffx <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/'
  tmp.replace <- '/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/deconvolution/genotype_based/demuxlet/'
  aggr.table[, demuxlet.tag:=str_replace_all(
    string=str_replace(string=molecule_h5, pattern=seq.pffx, replacement=tmp.replace),
    pattern='/count/|/outs/molecule_info.h5', replacement='/'
  )]
  aggr.table[, demuxlet.tag:=paste0(demuxlet.tag, basename(demuxlet.tag), '_', demuxlet.sffx, '/ClassificationPerCell.csv')]
  aggr.table[, hto.tag:=demuxlet.tag]
  tmp.check <- aggr.table[, all(file.exists(demuxlet.tag))]
  if(!tmp.check){
    tmp.data <- aggr.table[!file.exists(demuxlet.tag), demuxlet.tag]
    tmp.err <- paste0('Unexpected error. Next files do not exist:\n', paste0(tmp.data, collapse='\n'), '\n')
    stop(tmp.err)
  }
  # Make a copy of the current seurat object and modify seurat object's metadata to keep only base information.
  alt.obj <- seurat.obj
  cols.to.keep <- c('orig.ident', 'nCount_RNA', 'nFeature_RNA')
  alt.obj@meta.data <- alt.obj@meta.data[, cols.to.keep]
  # Annotate new seurat object.
  alt.obj <- annotate.seurat.obj(seurat.obj=alt.obj, aggr.table=aggr.table, lane.tag.id='lane.id.tag', chrom.batch.tag.id='chrom.batch.tag', seq.batch.tag.id='seq.batch.tag', hto.tag.id='demuxlet.tag', overlay.tag.id='donor.tag', mult.ann.per.lane=TRUE, exceded.batches=FALSE)
  alt.obj <- alt.obj@meta.data[, 'demuxlet.tag', drop=FALSE]
  # Merge previous donor data and new data.
  tmp.data <- merge(x=seurat.obj@meta.data, y=alt.obj, by='row.names', sort=FALSE, all=TRUE)
  row.names(tmp.data) <- tmp.data$Row.names; tmp.data$Row.names <- NULL
  # Modify classes provided by demuxlet to match those obtained based on the HTO data.
  demuxlet.vals <- tmp.data[, 'demuxlet.tag']
  demuxlet.vals[demuxlet.vals == 'DBL' | demuxlet.vals == 'AMB'] <- NA
  demuxlet.vals <- str_replace(string=demuxlet.vals, pattern='^DLCP[0]*', replacement='')
  # Then create the final/full donor ID
  tmp.data$full.donor.id.tag <- tmp.data$donor.id.tag
  tmp.cells <- is.na(tmp.data$donor.id.tag)
  tmp.data[tmp.cells, 'full.donor.id.tag'] <- demuxlet.vals[tmp.cells]
  # Replace meta data.
  tmp.check <- all(row.names(tmp.data) == Cells(seurat.obj))
  if(!tmp.check) stop('Unexpected error when creating the final donor ID tag.\n')
  seurat.obj@meta.data <- tmp.data
  return(seurat.obj)
}


# 11 ------------------------------------------------------------------->
# Name: Preprocess seurat object and provide meta data as well.
# Description:
# Given a 'raw' seurat object (i.e., as output by the basic seurat analysis pipeline), this function aims to preprocess it through next basic steps:
#   Definition of gene signatures.
#   Set factors properly.
#   Retrieve metadata.
#   Infer polyreactivity tag if needed.
# If you require to preprocess a given seurat object through steps other than the ones described above (i.e., specific for an object), do so before calling this function so that the metadata can include the appropriate features.
# Arguments ------------------------->
#
# Function:

preprocess.obj <- function(seurat.obj, set.lab, clusts.lab, is.ensembl=TRUE, donor.meta.data, do.tcr.data=FALSE, get.metadata=TRUE){
  # @ For eQTL analysis, rescue cells without HTO data.
  seurat.obj <- rep.hto.negatives.with.demuxlet(seurat.obj=seurat.obj, set.lab=set.lab, tissue.type=ifelse(test=str_detect(string=set.lab, pattern='tumor'), yes='Tumor-lung', no='Normal-lung'), aggr.pffx.lab='R24_Cancer_Batches-1-to-20_')
  # @ Further include the available donor demographics and clinical information.
  tmp.data <- seurat.obj@meta.data; tmp.data$cell <- Cells(seurat.obj)
  tmp.data <- merge(x=tmp.data, y=donor.meta.data, by.x='full.donor.id.tag', by.y='donor.id.tag', sort=FALSE, all.x=TRUE, all.y=FALSE)
  row.names(tmp.data) <- tmp.data$cell; tmp.data$cell <- NULL
  tmp.data <- tmp.data[Cells(seurat.obj), ]
  seurat.obj@meta.data <- tmp.data
  # @ Definition of signatures in seurat object.
  seurat.obj <- get.module.scores(seurat.obj=seurat.obj, module.feats.dir=module.sigs.path, signal.thold=0, is.ensembl=is.ensembl, reports.path=NULL)
  # @ Set factors appropriately
  # Clusters.
  new.lvls <- mixedsort(unique(as.character(seurat.obj@meta.data[, clusts.lab])))
  seurat.obj@meta.data[, clusts.lab] <- factor(x=as.character(seurat.obj@meta.data[, clusts.lab]), levels=new.lvls)
  # @ Take general metadata from the object.
  meta.data <- as.data.table(cbind(seurat.obj@meta.data, seurat.obj@reductions$umap@cell.embeddings))
  meta.data$cell <- Cells(seurat.obj)
  # @ Further tags.
  # Cell barcode.
  meta.data$cell <- as.character(meta.data$cell)
  # If required, update clone size tag solely by taking into acccount the overlap of both the gene expression and TCR datasets.
  if(do.tcr.data){
    tmp.data <- meta.data[!is.na(clonotype.tag), .(new.clon.size=.N), by=clonotype.tag]
    meta.data <- merge(x=meta.data, y=tmp.data, by='clonotype.tag', sort=FALSE, all=TRUE)
  }
  #   Sanity checks before proceeding to update actual seurat object.
  tmp.check.1 <- all(meta.data[, cell]==Cells(seurat.obj))
  tmp.check.2 <- if(do.tcr.data) meta.data[is.na(clonotype.tag), all(is.na(new.clon.size))] else TRUE
  tmp.check <- tmp.check.1 & tmp.check.2
  if(!tmp.check) stop('Unexpected error.\n')
  if(do.tcr.data) seurat.obj@meta.data[, 'clon.size.tag'] <- meta.data[, new.clon.size]
  # @ Output objects of interest.
  to.output <- if(get.metadata) list(seurat.obj, meta.data) else seurat.obj
  return(to.output)
}

# 12 ------------------------------------------------------------------->
# Name: Get publication figure.
# Description:
# Provided a ggplot, the name of a file to save it to and any further pdf function arguments, this function will output the publication quality pdf of such a plot in both versions: 1) 'Complete' version: version with all labels; and 2) 'Blank' version: version to be used in the actual paper figures file.
# Arguments ------------------------->
#
# Function:

publish.plot <- function(tmp.ggplot, output.path, file.name, type='pdf', blank.comp, do.legend=FALSE, stat.cor=FALSE, cor.group=NULL, repel.data=NULL, repel.label=NULL, do.rotate=FALSE, ...){
  # Blank version.
  tmp.file.name <- if(type=='pdf') paste0(output.path, '/', file.name, '.B.pdf') else paste0(output.path, '/', file.name, '.B.tiff')
  if(type=='pdf') pdf(file=tmp.file.name, ...) else tiff(file=tmp.file.name, height=1400, width=1400)
  print(tmp.ggplot + blank.comp)
  dev.off()
  # Complete version.
  tmp.ggplot <- if(stat.cor){
    if(!is.null(cor.group)) tmp.ggplot + stat_cor(aes_string(group=cor.group)) else tmp.ggplot + stat_cor()
  }else{
    to.check <- !is.null(repel.data) & !is.null(repel.label)
    if(to.check){
      tmp.ggplot + ggrepel::geom_text_repel(data=repel.data, aes_string(label=repel.label), color='black')
    }else{
      tmp.ggplot
    }
  }
  tmp.file.name <- if(type=='pdf') paste0(output.path, '/', file.name, '.C.pdf') else paste0(output.path, '/', file.name, '.C.tiff')
  tmp.ggplot <- if(do.rotate) tmp.ggplot + theme_bw() + theme(axis.text.x=element_text(angle=45)) else tmp.ggplot + theme_bw()
  if(type=='pdf') pdf(file=tmp.file.name, ...) else tiff(file=tmp.file.name, height=1400, width=1400)
  print(tmp.ggplot)
  dev.off()
  # Output legend independently.
  if(do.legend){
    tmp.ggplot <- get_legend(p=tmp.ggplot + theme(legend.position='right', legend.background=element_blank(), legend.key=element_blank(), legend.text=element_blank(), legend.title=element_blank()))
    tmp.file.name <- if(type=='pdf') paste0(output.path, '/', file.name, '.L.pdf') else paste0(output.path, '/', file.name, '.L.tiff')
    if(type=='pdf') pdf(file=tmp.file.name, ...) else tiff(file=tmp.file.name, height=1400, width=1400)
    print(as_ggplot(tmp.ggplot))
    dev.off()
  }
}

# 13 ------------------------------------------------------------------->
# Name: Plot proportions.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tags set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of one or multiple tags (filter.tags). See details below as to how to indicate filtring based on multiple tags.
# 5. Generate a ggplot depicting proportions of cells per group of a given metadata column, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided). The type of plot is defined through argument 'plot.type'.
# Arguments ------------------------->
#   seurat.obj - seurat object. The tags indicated in any of the arguments of this function must be already defined. No default.
#   groups.tag - Tag defined in the seurat onject's metadata. Defines the groups to depict cell fractions of.
#   na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
#   groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
#   groups.label - Character. When a non-NULL value is provided, indicates that the groups of interest provided throug argument 'groups.of.int' should be in turn merged into a single group and the value given through this argument will be used as a label for such a single group. Otherwise, a NULL value indicates that the groups of interest shouldn't be merged. Default: NULL.
#   filter.tags - Chracter vector, which should list the tag or tags in the seurat object metadata whose groups will be used to filter out cells. When filtering should be applied based on multiple tags, their groups to be considered to keep or exclude cells must be passed as chracter vectors embedded within a list; then, this character vector must have the same length as that of the list provided through the argument 'groups.to.filter'. No cell filtering is applied when it's set to NULL. Default: NULL.
#   groups.to.filter - List of vectors for the set of tag groups that'll be kept ot removed according to the argument 'keep' value. Must have the same length as the vector provided in 'filter.tags'. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
#   color.vals - Set of Groups' colors to be used as argument 'values' is 'scale_fill_manual'; thus, must be a character vector of color IDs with names corresponding to the groups of 'groups.tag'.
#   file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

plot.props <- function(
  # Input data and type of plot to output.
  seurat.obj, plot.type='fill',
  # Info to group and filter data.
  groups.tag, na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE,
  # Specificities for plotting.
  color.vals=NULL, file.name=NULL
){
  # ---> Arguments and checks.
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # ---> Fetch data.
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
    tmp.data <- cbind(data.table(groups.tag=seurat.obj@meta.data[, groups.tag]), tmp.data)
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='Rest']
  # Similarly, merge groups of interest if their label is provided.
    if(!is.null(groups.label)){
      tmp.data[groups.tag %in% groups.of.int, groups.tag:=as.character(groups.label)]
      # Set factors leaving 'Rest' label right at the end of the levels' order.
      new.lvls <- c(as.character(groups.label), 'Rest')
      tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
    }
  }
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Proportions plot.
  # Get ggplot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=0, fill=groups.tag)) + geom_bar(position='fill', alpha=0.8) + labs(x='', y='Proportion', fill='Group')
  # Adjust color scale.
  if(!is.null(color.vals)){
    tmp.check <- all(names(color.vals) %in% tmp.data[, unique(groups.tag)]) & all(tmp.data[, unique(groups.tag)] %in% names(color.vals))
    if(tmp.check){
      tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=color.vals)
    }else{
      tmp.warning <- 'Set of colors for the groups to be displayed were provided, but names in vector do not match with those to be displayed (from seurat object), so they will not be used.\n'
      warning(tmp.warning)
    }
  }
  if(plot.type=='donut') tmp.ggplot <- tmp.ggplot + coord_polar(theta = "y") + xlim(c(-2, 0.5)) + theme_void()
  if(plot.type=='pie') tmp.ggplot <- tmp.ggplot + coord_polar(theta = "y") + theme_void()
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}

# 14 ------------------------------------------------------------------->
# Name: Depict clonotype sharing.
# Description:
# Provided the metadata of a single-cell dataset (ideally, from a seurat object) and a tag defined in it, this function will depict the clonotype sharing among the values defined for the tag of interest. SHaring is depicted at two levels: 1. Clonotypes (venn diagram) and 2. Cells (bar plot depicting proportions).
# Arguments ------------------------->
#   meta.data - Single-cell dataset's metadata. May be a data table or not, but must be a dataframe.
#   tag.of.int - Tag defined in the seurat object's metadata. Sharing will be defined according to its groups.
#   size.thold -
#   reports.path - Absolute path.
#   file.preffix - File preffix. Appropriate extensions and suffix description will be appended to it for the output files.
#   fill.cols - Vector of colors for every group defined for the tag of interest in metadata.
#   sharing.cols - Colors for groups 'Shared' and 'Unique'

depict.clone.sharing <- function(
  meta.data, tag.of.int, size.thold=1, # Input data.
  filter.tag=NULL, groups.to.filter=NULL,
  reports.path, file.preffix, # About the output.
  fill.cols=NA, sharing.cols=c(Shared='#cccc00', Unique='#0074e6')
){ # Aesthetics.
  # @ Sharing according to clonotypes.
  # Harvest data.
  tmp.meta <- meta.data[!is.na(clonotype.tag)]
  if(!is.null(filter.tag)){
    if(!filter.tag %in% colnames(tmp.meta)) stop('Filter tag was provided but not found in the metadata.\n')
    tmp.meta[, tmp.filter.tag:=get(filter.tag)]
    tmp.check <- all(groups.to.filter %in% tmp.meta[, tmp.filter.tag])
    if(!tmp.check) stop('Not all of the groups indicated to filter forward were found in the metadata object.\n')
    tmp.meta <- tmp.meta[tmp.filter.tag %in% groups.to.filter]
    tmp.meta[, tmp.filter.tag:=NULL]
  }
  # Define data to depict.
  tmp.meta[, tag.of.int:=get(tag.of.int)]
  tmp.vals <- tmp.meta[, unique(tag.of.int)]
  tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.val){
    to.return <- tmp.meta[tag.of.int==tmp.val, .(count=.N), by=.(clonotype=clonotype.tag)][count>size.thold, clonotype] # This may be changed to adjust for the minimum size for a clonotype to be taken into account.
    return(to.return)
  })
  names(tmp.data) <- tmp.vals
  # Plot and output.
  tmp.file.name <- paste0(reports.path, '/', file.preffix, '-Clonotypes.C.tiff')
  venn.diagram(x=tmp.data, filename=tmp.file.name, imagetype='tiff', fill=fill.cols, col=fill.cols)
  tmp.file.name <- paste0(reports.path, '/', file.preffix, '-Clonotypes.B.tiff')
  venn.diagram(x=tmp.data, filename=tmp.file.name, imagetype='tiff', fill=fill.cols, col=fill.cols, cex=0, cat.cex=0, ext.text=FALSE)
  # Remove temporary files.
  tmp.files <- list.files(path=reports.path, pattern='\\d+-\\d+-\\d+\\.log$', full.names=TRUE)
  file.remove(tmp.files)
  # @ Sharing according to cells.
  # Harvest data.
  tmp.data <- Reduce(x=tmp.data, f=intersect)
  tmp.meta[, tcr.sharing.tag:='Unique']
  tmp.meta[clonotype.tag %chin% tmp.data, tcr.sharing.tag:='Shared']
  # For caption describing number of cells per sharing group per tag group.
  tmp.data.2 <- tmp.meta[, .(cell.count=.N), by=.(group=tag.of.int, sharing.status=tcr.sharing.tag)]
  tmp.file.name <- paste0(reports.path, '/', file.preffix, '-Cells.csv')
  fwrite(file=tmp.file.name, x=tmp.data.2, na=NA, quote=FALSE)
  # Get plot.
  tmp.ggplot <- ggplot(data=tmp.meta, aes(x=tag.of.int, fill=tcr.sharing.tag)) +
    geom_bar(position='fill', width=0.6) +
    scale_y_continuous(expand=expansion(add=c(0, 0))) +
    scale_fill_manual(values=sharing.cols) +
    labs(x='', y='Cell fraction', fill='')
  # Output
  publish.plot(tmp.ggplot=tmp.ggplot, output.path=reports.path, file.name=paste0(file.preffix, '-Cells'), blank.comp=blank.complement.1, do.legend=TRUE)
}

# 15 ------------------------------------------------------------------->
# Name: Get metrics for ranking features.
# Description:
# Given a seurat object and a tag listing a value of interest (defined in the object's metadata and provided to this function as 'tag.of.int'), this function will calculate the metrics to rank the object's features taking the value of interest (specified in argument 'val.of.int') as the reference group and the rest of the cells as the alternative group. Ranking metrics that can be calculated by this function follow (calling name accepted by the function right after each of them):
#   Signal to noise ratio: 'signal.to.noise'
#   T test stat: 't.test'
#   Ratio of classes' means: 'ratio.of.classes'
#   Difference of classes' means: 'diff.of.classes'
#   Log2 ratio of classes' means: 'log.ratio.of.classes'
# For further information about such metrics, visit: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
# Arguments ------------------------->
# seurat.obj - A seurat object with features for which to calculate metrics.
# tag.of.int - Tag name (as defined in the seurat object's metadata) listing value of interest (to be taken as reference group).
# val.of.int - Value of interest whose cell'll define the group of interets.
# metric - Metric for ranking features to be calculated. Accepted metrics are: signal.to.noise, t.test, ratio.of.classes, diff.of.classes, log.ratio.of.classes
# Function:

get.rank.metrics <- function(seurat.obj, tag.of.int, val.of.int, metric='signal.to.noise'){
  # Preflights.
  tmp.check <- tag.of.int %in% colnames(seurat.obj@meta.data)
  if(!tmp.check){ tmp.error <- paste0('When getting ranking metrics for tag value ', val.of.int, ' as reference, we could not find tag ', tag.of.int, '. Make sure to provide a valid -already defined in seurat object- tag of interest.\n'); stop(tmp.error) }
  tmp.check <- val.of.int %in% seurat.obj@meta.data[, tag.of.int]
  if(!tmp.check){ tmp.error <- paste0('Value ', val.of.int, ' could not be found at all as defined for tag ', tag.of.int, '. Make sure to provide valid tag values for such a tag.\n'); stop(tmp.error) }
  accepted.metrics <- c('signal.to.choice', 't.test', 'ratio.of.classes', 'diff.of.classes', 'log.ratio.of.classes')
  tmp.check <- as.character(metric) %chin% accepted.metrics
  if(tmp.check){ tmp.error <- paste0('Metric requested, ', metric, ', is not a valid one. Metrics accepted by this fucntion are: ', paste0(accepted.metrics, collapse=', '), '. Make you have typed your option correctly.\n'); stop(tmp.error) }
  # Group cells according to group of interest.
  ref.group <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.int]==val.of.int & !is.na(seurat.obj@meta.data[, tag.of.int])]
  alt.group <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.int]!=val.of.int | is.na(seurat.obj@meta.data[, tag.of.int])]
  # Kill process by returning an NA value if reference group size is too large or too small.
  min.size <- min(length(ref.group), length(alt.group))
  tmp.check <- min.size < 3
  if(tmp.check){ tmp.warning <- paste0('Process was not finished because cells tagged with value ', val.of.int, ' compose a group with a size that is too different compared to the size of the group composed by the rest of the cells. Returning NA, instead.\n'); warning(tmp.warning); return(NA) }
  # Calculate genes' basic stats.
  # @ Count of samples.
  ref.n <- length(ref.group)
  alt.n <- length(alt.group)
  # @ Mean
  ref.means <- Matrix::rowMeans(GetAssayData(object=seurat.obj, slot='data', assay='RNA')[, ref.group])
  alt.means <- Matrix::rowMeans(GetAssayData(object=seurat.obj, slot='data', assay='RNA')[, alt.group])
  # @ Standard deviation.
  ref.sds <- apply(X=GetAssayData(object=seurat.obj, slot='data', assay='RNA')[, ref.group], MARGIN=1, FUN=sd)
  alt.sds <- apply(X=GetAssayData(object=seurat.obj, slot='data', assay='RNA')[, alt.group], MARGIN=1, FUN=sd)
  # @ Set a minimum for standard deviation (per GSEA general recommendations -https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html-) and substitute all NaN values with such a minimum value.
  tmp.means <- ref.means; tmp.means[tmp.means==0] <- 1
  ref.min.sd <- 0.2 * abs(tmp.means)
  ref.sds[ref.sds<ref.min.sd] <- ref.min.sd[ref.sds<ref.min.sd]
  tmp.means <- alt.means; tmp.means[tmp.means==0] <- 1
  alt.min.sd <- 0.2 * abs(tmp.means)
  alt.sds[alt.sds<alt.min.sd] <- alt.min.sd[alt.sds<alt.min.sd]
  # Get metric according to the user's choice.
  to.output <- switch(
    EXPR=metric,
    signal.to.noise={
      cat('Calculating signal to noise ratio...\n')
      (ref.means - alt.means) / (ref.sds + alt.sds)
    },
    t.test={
      cat('Calculating T test...\n')
      denominator <- sqrt(((ref.sds ** 2) / ref.n) + ((alt.sds ** 2) / alt.n))
      (ref.means - alt.means) / denominator
    },
    ratio.of.classes={
      cat('Calculating ratio of classes...\n')
      ref.means / alt.means
    },
    diff.of.classes={
      cat('Calculating differences of classes...\n')
      ref.means - alt.means
    },
    log.ratio.of.classes={
      cat('Calculating logged ratio of classes...\n')
      log2(x=(ref.means / alt.means))
    }
  )
  return(to.output)
}

# 16 ------------------------------------------------------------------->
# Name: Do FGSEA.
# Description:
# Given a seurat object and a tag listing a value of interest (defined in the object's metadata and provided to this function as 'tag.of.int'), this function will apply FGSEA comparing the signature enrichment between a set of values defined for the tag of interest (for each one independently) and the rest of the cells. Ranking metrics are calculated independently (see function above, 'get.rank.metrics') and same values are valid for this function.
# NOTE: If you wish to use this function for the development of other papers, pay attention to how files are being output here (see argument 'files.preffix').
# Arguments ------------------------->
# metrics.src - A seurat object with features for which to calculate metrics. In its metadata the groups of interest are defined under the name provided through argument 'tag.of.int'.
# tag.of.int - Tag name (as defined in the seurat object's metadata) listing value of interest (to be taken as reference group).
# vals.to.depict - Character, lists the sets of values defined for the tag of interest for which an output file (GSEA plot) should be output.
# metric - Metric for ranking features to be calculated. Accepted metrics are: signal.to.noise, t.test, ratio.of.classes, diff.of.classes, log.ratio.of.classes (see function above, 'get.rank.metrics')
# files.preffix - Preffix to add to every plot to be output by this function.
# Function:

do.fgsea <- function(metrics.src, modules.feats, tag.of.int, metric='signal.to.noise', output.path, vals.to.depict=NULL, do.heatmap=FALSE, full.reports=FALSE, files.preffix, signature.order=NULL, ref.order=NULL, col.metadata=NULL, ann.colors=ann.colors, just.get.res=FALSE){
  # ---> Predefinitions.
  files.preffix <- paste0(files.preffix, '_FGSEA-')
  # ---> Metrics calculation for expresion values.
  # Define tag values and obtain metrics for each of them taken as reference.
  tag.values <- mixedsort(unique(metrics.src@meta.data[, tag.of.int]))
  tag.values <- tag.values[!is.na(tag.values)]
  if(is.factor(tag.values)) tag.values <- as.character(tag.values)
  val.metrics <- lapply(X=tag.values, FUN=function(val.of.int) get.rank.metrics(seurat.obj=metrics.src, tag.of.int=tag.of.int, val.of.int=val.of.int, metric=metric))
  names(val.metrics) <- tag.values
  # Confirm values to be depicted that were provided are indeed defined within the seurat object.
  if(is.null(vals.to.depict) | do.heatmap){
    vals.to.depict <- tag.values
  }else{
    to.check <- all(vals.to.depict %in% tag.values)
    if(!to.check){
      to.check <- vals.to.depict %in% tag.values
      to.check <- vals.to.depict[!to.check]
      tmp.err <- paste0('Next values were requested to be depicted, but they have not been previously defined for tag \'', tag.of.int, '\':\n',
     paste0(to.check, collapse='\n'), '\nPlease make sure to provide valid values for the tag of interest.\n')
     stop(tmp.err)
    }
  }
  # Remove values for which metrics calculation was not successful.
  to.check <- is.na(val.metrics)
  to.check <- names(to.check)[to.check]
  if(length(to.check)>0){
    tmp.warning <- paste('Metrics calculation was not successful for next values because their group size was too small:\n\t', paste0(to.check, collapse='\t\n'), '\nThey will be removed from further processes.\n')
    warning(tmp.warning)
    val.metrics <- val.metrics[!is.na(val.metrics)]
    to.check <- length(val.metrics)<1
    if(to.check){ tmp.error <- paste0('Seems that metrics calculation failed for all of the possible value from the tag of interest (', tag.of.int, '), so the program cannot continue any further.\n'); stop(tmp.error)}
  }
  # Save calculated metrics.
  if(full.reports){
    to.output <- rbindlist(l=lapply(X=val.metrics, FUN=function(x) data.table(gene=names(x), metric=x)), use.names=TRUE, idcol='tag.value')
    to.output <- spread(data=to.output, key=tag.value, value=metric)
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'RankingMetricsPerTagValueAcrossGenes_Tag-', tag.of.int, '.csv')
    fwrite(file=tmp.file.name, x=to.output, na='NaN', quote=FALSE)
  }
  # ---> FGSEA
  cols.to.keep <- c('pathway', 'pval', 'padj', 'ES', 'NES', 'nMoreExtreme', 'size')
  fgsea.results <- lapply(X=vals.to.depict, FUN=function(tag.value){
    # FGSEA
    set.of.metrics <- val.metrics[[tag.value]]
    set.of.metrics <- set.of.metrics[!is.nan(set.of.metrics)]
    fgsea.results <- fgsea(pathways=modules.feats, stats=set.of.metrics, minSize=10, maxSize=500, nproc=total.workers, nperm=10000)
    # Save plain results.
    fgsea.results <- fgsea.results[, ..cols.to.keep]
    if(full.reports){
      tmp.file.name <- paste0(output.path, '/', files.preffix, 'FGSEAResultsForTagValue_Tag-', tag.of.int, '_Val-', tag.value, '.csv')
      fwrite(file=tmp.file.name, x=fgsea.results, na=NA, quote=FALSE)
    }
    # Return general results.
    return(fgsea.results)
  })
  names(fgsea.results) <- vals.to.depict
  # ---> Just get results
  if(just.get.res) return(fgsea.results)
  # ---> Individual enrichment plots.
  if(do.heatmap){
    # Tidy data.
    tmp.data <- rbindlist(l=fgsea.results, use.names=TRUE, idcol='ref.tag')
    tmp.data <- tmp.data[, .(ref.tag, signature=pathway, p.adj=padj, NES)]
    tmp.data[p.adj>0.05, NES:=NA]
    tmp.data[, p.adj:=NULL]
    tmp.data <- spread(data=tmp.data, key=ref.tag, value=NES)
    tmp.data <- as.data.frame(tmp.data, stringsAsFactors=FALSE)
    row.names(tmp.data) <- tmp.data$signature; tmp.data$signature <- NULL
    # Set color scale and breaks for heatmap.
    col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=100)
    mid.point <- which.min(abs(col.breaks - 0))
    # hmap.col.scale.1 <- colorRampPalette(c('mediumblue', 'lightblue', 'white'))(mid.point)
    hmap.col.scale.1 <- colorRampPalette(c('#233773', '#3553ae', '#72bcd4', '#ffffff'))(mid.point)
    # hmap.col.scale.2 <- colorRampPalette(c('white', 'yellow', 'orange', 'red'))(100-(mid.point+1))
    hmap.col.scale.2 <- colorRampPalette(c('#ffffff', '#ffff00', '#ffa500', '#ff664d', '#ff2500', '#b31a00'))(100-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Set order of data.
    if(!is.null(signature.order)){
      if(is.null(names(signature.order))) names(signature.order) <- signature.order
      tmp.data <- tmp.data[signature.order, ]
      signature.order <- setNames(object=names(signature.order), nm=signature.order)
      row.names(tmp.data) <- signature.order[row.names(tmp.data)]
    }
    if(!is.null(ref.order)){
      if(is.null(names(ref.order))) names(ref.order) <- ref.order
      tmp.data <- tmp.data[, ref.order]
      ref.order <- setNames(object=names(ref.order), nm=ref.order)
      colnames(tmp.data) <- ref.order[colnames(tmp.data)]
    }
    # @ Heatmap.
    # Plot 'complete' version
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.C.pdf')
    pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, show_rownames=TRUE, show_colnames=TRUE, filename=tmp.file.name, na_col='#b3b3b3')
    # Plot 'blank' version
    tmp.file.name <- paste0(output.path, '/', files.preffix, 'SignaturesNESHeatmap.B.pdf')
    pheatmap(mat=tmp.data, color=hmap.col.scale, breaks=col.breaks, scale='none', cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=col.metadata, annotation_colors=ann.colors, legend=TRUE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE, show_rownames=FALSE, show_colnames=FALSE, filename=tmp.file.name)
  }else{
    for(tag.value in vals.to.depict){
      for(tmp.module in names(modules.feats)){
        # Check pathway was actually processed.
        tmp.check <- tmp.module %in% fgsea.results[[tag.value]][, pathway]
        if(!tmp.check) next
        # Retrieve stats, module signatures and further labels.
        tmp.lab <- toupper(str_replace_all(string=tmp.module, pattern='_', replacement=' '))
        tmp.pathway <- modules.feats[[tmp.module]]
        tmp.stats <- val.metrics[[tag.value]]
        tmp.stats <- tmp.stats[!is.nan(tmp.stats)]
        fgsea.stats <- fgsea.results[[tag.value]][pathway==tmp.module, round(x=c(ES, NES, padj), digits=6)]
        tmp.caption <- paste0('ES: ', fgsea.stats[1], '; NES: ', fgsea.stats[2], '; adj. P value: ', fgsea.stats[3])
        # Get plot.
        tmp.ggplot <- plotEnrichment(pathway=tmp.pathway, stats=tmp.stats) + labs(title=tmp.lab, x='Rank', y='Enrichment score', caption=tmp.caption)
        # Output enrichment plot.
        tmp.file.name <- paste0(files.preffix, 'EnrichmentPlotForSignature-', tmp.module, '_Tag-', tag.of.int, '_Val-', tag.value)
        publish.plot(tmp.ggplot=tmp.ggplot, output.path=output.path, file.name=tmp.file.name, type='pdf', blank.comp=blank.complement.1, do.legend=FALSE, height=4, width=4)
      }
    }
  }
}

# 17 ------------------------------------------------------------------->
# Name: read.seurat.input.data
# Description:
# Given a file/path name, this function figures out what its format is to work as a data object appropriate to become a seurat object. It looks for three formats
# a) TSV with extensions .tsv and .txt
# b) CSV with extension csv
# c) H5 with extension h5
# The default is for it to be taken as a path to cellranger output where matrix, features and cells data is saved.
# The function returns the file read and saved to an object.
# Arguments ------------------------->
# file - Character value to be assessed for any of the different formats described.
# feature.id - Character with possible values 'name' or 'ensembl'. Used just in case a cellranger output is provided (either from features.csv file or from h5 matrix). Default, name.

read.10X.data <- function(file, feature.id='name'){
  if(!is.null(dim(file))) stop('file argument must be a character vector storing in its first index the name of the file to be read.')
  file <- file[1]
  # ---> Check feature.id option value.
  # Valid values: 'name', 'ensembl'
  if(!(feature.id=='name'|feature.id=='ensembl')) stop('Non valid entry for feature.id argument.')
  # ---> Read data.
  #   ---> Normal table
  # Tab-separated values.
  if(grepl(x=file, pattern='\\.txt$', ignore.case=TRUE, perl=TRUE) | grepl(x=file, pattern='\\.tsv', ignore.case=TRUE, perl=TRUE)) return(read.table(file=file))
  # Comma-separated values.
  if(grepl(x=file, pattern='\\.csv$', ignore.case=TRUE, perl=TRUE)) return(read.csv(file=file, row.names=1))
  #   ---> Cellranger output.
  # Case h5 format.
  if(grepl(x=file, pattern='.h5$', ignore.case=TRUE, perl=TRUE)){
    feature.id <- feature.id=='name'
    return(Seurat::Read10X_h5(filename=file, use.name=feature.id))
  }
  # Directory with 10X output (default).
  return(Seurat::Read10X(data.dir=file, gene.column=ifelse(test=feature.id=='name', yes=2, no=1)))
}

# 18 ------------------------------------------------------------------->
# Name: Raw seurat object.
# Description:
# Provided the ID of a seurat object, this function will provide its raw seurat object with annotations (i.e., previous to any filtering step).
# Arguments ------------------------->
# obj.of.int - Seurat object ID of interest.
# Function:

get.raw.obj <- function(obj.of.int){
  # @ File definitions and others.
  data.file <- paste0(gen.data.path, '/RawData_', obj.of.int)
  aggr.table.file <- paste0(gen.data.path, '/ForAnns_', obj.of.int, '.csv')
  if(!file.exists(aggr.table.file)) stop('Aggregation table file does not exist.\n')
  donors.data.file <- paste0(gen.data.path, '/DonorsMetadata_', obj.of.int, '.csv')
  donors.merge.by <- 'donor.tag'
  # tags.criteria.file <- paste0(gen.data.path, '/Filt_', obj.of.int, '.csv')
  # @ Load data.
  # Raw data (create seurat object).
  data <- read.10X.data(file=data.file, feature.id='name') # Create seurat object.
  seurat.obj <- CreateSeuratObject(counts=data, project=obj.of.int, min.cells=0)
  rm(data)
  # Annotations table.
  aggr.table <- read.csv(file=aggr.table.file, stringsAsFactors=FALSE)
  donors.data <- read.csv(file=donors.data.file, stringsAsFactors=FALSE)
  # Subset criteria.
  # if(!is.null(tags.criteria.file)) tags.criteria <- read.csv(file=tags.criteria.file, stringsAsFactors=FALSE, row.names=1) else tags.criteria <- NULL
  # @ Annotations
  seurat.obj <- annotate.seurat.obj(seurat.obj=seurat.obj, aggr.table=aggr.table, lane.tag.id=lane.tags[[obj.of.int]], chrom.batch.tag.id='chrom.batch.tag', seq.batch.tag.id='seq.batch.tag', hto.tag.id='hashtag.tag', overlay.tag.id='donor.tag', mult.ann.per.lane=TRUE)
  # @ Annotations based on donors' data.
  if(!(donors.merge.by %in% colnames(donors.data) & donors.merge.by %in% colnames(seurat.obj@meta.data))) stop('Donors metadata was provided but not an appropriate value to add further annotations (i.e., value was not part of seurat annotations.)')
  # ---> Find per-cell further annotations.
  tmp.df <- data.frame(col.1=seurat.obj@meta.data[, donors.merge.by], barcode=Cells(seurat.obj), stringsAsFactors=FALSE)
  colnames(tmp.df)[1] <- donors.merge.by
  tmp.merge <- merge(x=tmp.df, y=donors.data, by=donors.merge.by, all.x=TRUE)
  rownames(tmp.merge) <- tmp.merge$barcode
  tmp.merge$barcode <- NULL
  tmp.merge <- tmp.merge[Cells(seurat.obj), ]
  if(!all(rownames(tmp.merge) == Cells(seurat.obj))) stop('Something went wrong while adding patients\' annotations.')
  if(is.null(dim(tmp.merge[, colnames(tmp.merge)!=donors.merge.by]))){
    tmp.col.name <- colnames(tmp.merge)[colnames(tmp.merge)!=donors.merge.by]
    tmp.merge <- data.frame(tmp.merge[, colnames(tmp.merge)!=donors.merge.by])
    colnames(tmp.merge) <- tmp.col.name
  }else{
    tmp.merge <- tmp.merge[, colnames(tmp.merge)!=donors.merge.by]
  }
  tmp <- cbind(seurat.obj@meta.data, tmp.merge)
  # ---> Add all annotations.
  seurat.obj@meta.data <- tmp
  rm(tmp)
  # @ Return object.
  return(seurat.obj)
}

# 19 ------------------------------------------------------------------->
# Name: Annotate seurat object
# Dependencies:
#   Seurat and stringr
# Functions dependencies:
#   None.
# Description:
#    Given an aggregation table with annotations regarding each library (i.e., this is intended for seurat object created based on aggregations), this function will provide the appropriate annotations for each cell in the metadata section. Possible annotations are as follows:
#   HTO tags:
#     Description: Hashtag class per cell based on a file for each independent library, as reported by any method applied over hashtag data.
#     Value in aggr table: Absolute path to file describing class per cell.
#   HTO tags:
#     Description: Hashtag class per cell based on a file for each independent library.
#     Value in aggr table: Absolute path to file describing class per cell.
#   Chromium batch tag:
#     Description: Library ID.
#     Value in aggr table: Library ID.
#   Lane tag:
#     Description: Other info conveyed for the library ID. Could be a single value or different tags separated by semicolon.
#     Value in aggr table: Tag value or values separated by semicolon.
#   Chromium batch tag:
#     Description: Library ID.
#     Value in aggr table: Library ID.
#   Sequencing batch tag:
#     Description: Sequencing batch.
#     Value in aggr table: Sequencing batch.
#   Overlay tag:
#     Description: Value conveyed based on the combination of both tags, chromium batch and lane tag.
#     Value in aggr table: None. Added in case the base tags are provided.
# Arguments ------------------------->
# All argument values are taken from the main program.
#   seurat.obj, norm.method (which should be LogNormalize), FVFs.method, mean.cutoff, prop.cutoff, feats.for.dsa.
# Value ----------------------------->
# List containing:
#   Idx. 1) Processed seurat object.
#   Idx. 2) Chracter vector with final variable features.

annotate.seurat.obj <- function(seurat.obj, aggr.table, lane.tag.id, chrom.batch.tag.id, seq.batch.tag.id, hto.tag.id, overlay.tag.id, mult.ann.per.lane, exceded.batches=FALSE){
  ### ---------------------- Data preprocessing ----------------------- ###
  # ---> Samples suffix per cell.
  # Get sample-specific suffix of cell IDs.
  cell.sffxs <- as.integer(str_extract(string=Cells(seurat.obj), pattern='\\d+$'))
  # ---> Batches post-aggregation.
  # Find amount of 10x batches aggregated based on cell barcode suffixes and make sure it's the same amount of batch definitions provided.
  tenx.batches <- as.integer(unique(str_extract(string=Cells(seurat.obj), pattern='\\d+$')))
  if(length(tenx.batches)!=nrow(aggr.table) & !exceded.batches) stop('Not same amount of batches identified in Seurat object as the ones provided in table.')
  # ---> ClassificationPerCell files from aggregation table.
  # Just in case there's HTO data provided.
  if('hto.tag' %in% colnames(aggr.table)){
    hto.tag.vals.list <- lapply(X=aggr.table$hto.tag, FUN=function(hto.tag.file){
      if(is.na(hto.tag.file)){
        hto.tag.vals <- NULL
      }else{
        hto.tag.vals <- read.csv(file=hto.tag.file, stringsAsFactors=FALSE)
        colnames(hto.tag.vals) <- c('barcode', 'hashtag')
      }
      return(hto.tag.vals)
    })
    # Names designated as batch order.
    names(hto.tag.vals.list) <- as.character(1:length(hto.tag.vals.list))
    # Remove NULL values (batches with no HTO info.)
    hto.tag.vals.list <- hto.tag.vals.list[!sapply(X=hto.tag.vals.list, FUN=is.null)]
    if(length(hto.tag.vals.list)==0){
      hto.tag.vals <- NULL
    }else{
      # Unlist sets of HTO tag values and assing identifier according to batch. Then, names will be as in seurat object post-aggregation, composed by barcode and batch.
      # First, get complete cell IDs.
      cell.ids <- unlist(lapply(X=names(hto.tag.vals.list), FUN=function(batch.sfx){
        cell.id <- str_extract(string=hto.tag.vals.list[[batch.sfx]]$barcode, pattern='[ACTG]+')
        cell.id <- paste0(cell.id, '-',batch.sfx)
        return(cell.id)
      }))
      # Then get hashtag vaues unlisted and assign names.
      hto.tag.vals <- unlist(lapply(X=hto.tag.vals.list, FUN=function(hto.tag.vals) return(hto.tag.vals$hashtag)))
      names(hto.tag.vals) <- cell.ids
    }
  }
  ### ---------------------------- Lane tag --------------------------- ###
  if('lane.tag' %in% colnames(aggr.table)){
    # If there are multiple annotations.
    if(mult.ann.per.lane){
      # Tag values.
      lane.tags <- as.data.frame(str_split(string=aggr.table$lane.tag, pattern=';', simplify=TRUE))
      # Correct for NA vals.
      lane.tags[lane.tags=='NA'] <- NA
      # Tag IDs.
      lane.tag.id <- str_split(string=lane.tag.id, pattern=';', simplify=TRUE)[1, ]
      # Add one if it is not properly defined (best guess).
      if(length(lane.tag.id)!=ncol(lane.tags)){
        lane.tag.id <- paste(lane.tag.id[1], 1:ncol(lane.tags), sep='.')
        warning('No appropriate lane tag IDs were defined according to the number of values in aggr table. Next, number of values and number of IDs: ', paste(ncol(lane.tags), length(lane.tag.id), sep='/'))
      }
      # Then, set format.
      colnames(lane.tags) <- lane.tag.id
      lane.tag.vals <- lapply(X=1:nrow(lane.tags), FUN=function(tag.val) return(lane.tags[tag.val, ]))
      # Find lane tag value per cell.
      lane.val.per.cell <- rbindlist(lane.tag.vals[cell.sffxs])
      # Add lane tag value per cell.
      seurat.obj@meta.data <- cbind(seurat.obj@meta.data, lane.val.per.cell)
    }else{
    # If there's a single value.
      lane.tag.vals <- lapply(X=aggr.table$lane.tag, FUN=function(tag.val) return(tag.val))
      # Find lane tag value per cell.
      lane.val.per.cell <- unlist(lane.tag.vals[cell.sffxs])
      # Add lane tag value per cell.
      seurat.obj[[lane.tag.id]] <- lane.val.per.cell
    }
  }
  ### ---------------------------- HTO tag ---------------------------- ###
  if('hto.tag' %in% colnames(aggr.table)){
    if(is.null(hto.tag.vals)){
      seurat.obj[[hto.tag.id]] <- NA
    }else{
      classes.df <- data.frame(barcode.id=names(hto.tag.vals), class=hto.tag.vals, stringsAsFactors=FALSE)
      cells.df <- data.frame(barcode.id=Cells(seurat.obj), stringsAsFactors=FALSE)
      hto.tag.vals.per.cell <- merge(x=cells.df, y=classes.df, by='barcode.id', all.y=FALSE, all.x=TRUE, sort=FALSE)
      rownames(hto.tag.vals.per.cell) <- hto.tag.vals.per.cell$barcode.id
      hto.tag.vals.per.cell <- hto.tag.vals.per.cell[cells.df$barcode.id, ]
      # Add HTO tag value per cell.
      seurat.obj[[hto.tag.id]] <- hto.tag.vals.per.cell$class
    }
  }
  ### ----------------------- Chromium batch tag ---------------------- ###
  if('chrom.batch.tag' %in% colnames(aggr.table)){
    # Find batch tag value per cell.
    chrom.batch.tag.vals <- lapply(X=aggr.table$chrom.batch.tag, FUN=function(tag.val) return(tag.val))
    chrom.batch.val.per.cell <- unlist(chrom.batch.tag.vals[cell.sffxs])
    # Add batch tag value per cell.
    seurat.obj[[chrom.batch.tag.id]] <- chrom.batch.val.per.cell
  }
  ### ---------------------- Sequencing batch tag --------------------- ###
  if('seq.batch.tag' %in% colnames(aggr.table)){
    # Find batch tag value per cell.
    seq.batch.tag.vals <- lapply(X=aggr.table$seq.batch.tag, FUN=function(tag.val) return(tag.val))
    seq.batch.val.per.cell <- unlist(seq.batch.tag.vals[cell.sffxs])
    # Add batch tag value per cell.
    seurat.obj[[seq.batch.tag.id]] <- seq.batch.val.per.cell
  }
  ### ------------------------- Overlay tag --------------------------- ###
  # Tag provided by the combination of both, HTO and chromium batch tags.
  if(all(c('chrom.batch.tag', 'hto.tag') %in% colnames(aggr.table))){
    overlay.vals.per.cell <- paste(seurat.obj@meta.data[, hto.tag.id], seurat.obj@meta.data[, chrom.batch.tag.id], sep='-')
    # Find NA values from any independent column and report the same for overlay.
    overlay.vals.per.cell[is.na(seurat.obj@meta.data[, hto.tag.id]) | is.na(seurat.obj@meta.data[, chrom.batch.tag.id])] <- NA
    # Add overlay tag per cell.
    seurat.obj[[overlay.tag.id]] <- overlay.vals.per.cell
  }
  # ---> Return
  return(seurat.obj)
}

# 20 ------------------------------------------------------------------->
# Name: Merge QC results.
# Description:
# Given a set of metrics summary QC files, these will be merged depicting each sample as rownames.
# Arguments ------------------------->
# files.set - Character vector of absolute paths to qc summary files.
# Function:

merge.qcs <- function(files.set){
  # Get all tables to a list.
  qc.table.list <- lapply(X=files.set, FUN=function(tmp.file){
    tmp.qc.table <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
    samples.name <- basename(path=tmp.file)
    samples.name <- str_replace(string=samples.name, pattern="_metrics_summary\\.csv", replacement="")
    tmp.qc.table$Sample.ID <- samples.name
    return(tmp.qc.table)
  })
  # Merge them all together.
  qc.table.merge <- Reduce(x=qc.table.list, f=function(table.1, table.2) return(merge(table.1, table.2, all=TRUE)))
  # Add smaple names as row names.
  rownames(qc.table.merge) <- qc.table.merge$Sample.ID
  qc.table.merge$Sample.ID <- NULL
  # Then, return.
  return(qc.table.merge)
}

# 21 ------------------------------------------------------------------->
# Name: Extract aggr QCs.
# Description:
# Given a list obtained by reading the a metrics summary json file output by cellranger aggr (using fromJSON), this function will provide a the relevant information from such data to be reportted in the publication tables.
# Arguments ------------------------->
# json.data - List obtained by reading the a metrics summary json file output by cellranger aggr (using fromJSON)
# Value ----------------------------->
# List containing:
#   Idx. 1) 'gen.qc', data table with general metric values.
#   Idx. 2) 'libs.qc', data table listing the library-specific metric values.
# Function:

process.aggr.json <- function(obj.of.int, lib.name.pttn){
  # General definitions.
  json.data <- aggr.qcs[[obj.of.int]]
  gen.aggr.metrics <- c(
    post_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc = 'post.mean',
    post_normalization_total_reads = 'post.total',
    pre_normalization_multi_transcriptome_total_raw_reads_per_filtered_bc = 'pre.mean',
    pre_normalization_total_reads = 'pre.total'
  )
  aggr.metrics <- c(
    frac_reads_kept='frac.reads.kept',
    pre_normalization_cmb_reads_per_filtered_bc='cmb.reads.per.bc',
    pre_normalization_raw_reads_per_filtered_bc='raw.reads.per.bc'
  )
  aggr.metrics.names <- c(
    frac.reads.kept='Fraction of reads kept post-aggr',
    cmb.reads.per.bc='Average of reads confidently mapped to the genome per cell',
    raw.reads.per.bc='Average of raw reads per cell'
  )
  extract.metric.pttn <- paste0(names(aggr.metrics), collapse='|')
  rep.metric.pttn <- paste0('_', names(aggr.metrics), collapse='|')
  # Messy format as a data table.
  json.data <- unlist(json.data)
  json.data <- data.table(group=names(json.data), value=json.data)
  # Extract main metrics.
  gen.qc <- json.data[group %in% names(gen.aggr.metrics)]
  gen.qc[, group:=gen.aggr.metrics[group]]
  gen.qc[, value:=as.numeric(value)]
  # Extract metrics per library..
  json.data[, aggr.metric:=str_extract(string=group, pattern=extract.metric.pttn)]
  json.data[, sample.id:=str_replace(string=group, pattern=rep.metric.pttn, replacement='')]
  json.data[, sample.id:=str_replace(string=sample.id, pattern=paste0(obj.of.int, '\\.'), replacement='')]
  json.data <- json.data[!is.na(aggr.metric)]
  json.data <- json.data[str_detect(string=sample.id, pattern=lib.name.pttn)]
  json.data[, group:=NULL]
  json.data[, aggr.metric:=aggr.metrics[aggr.metric]]
  json.data[, value:=as.numeric(value)]
  json.data <- spread(data=json.data, key=aggr.metric, value=value)
  json.data[, cmb.reads.per.bc:=as.numeric(cmb.reads.per.bc)]
  json.data[, frac.reads.kept:=as.numeric(frac.reads.kept)]
  json.data[, raw.reads.per.bc:=as.numeric(raw.reads.per.bc)]
  return(list(gen.qc=gen.qc, libs.qc=json.data))
}

# 22 ------------------------------------------------------------------->
# Name: Get vdj gene usage information.
# Description:
# Provided the ID of an object, this function will load its 'filtered_clonotypes' file (saved to the general data folder with the pattern 'FilteredContigAnns_<OBJECT OF INTEREST>.csv') and provide the vdj gene usage information of each clonotype (i.e., the genes TRBV, TRBJ, TRAV and TRAJ for each clonotype that are best supported according to the amount of contigs annotated with a given gene). When multiple genes (of any kind) are equally well supported, all genes are provided.
# Arguments ------------------------->
# obj.of.int - ID of the object of interest.
# Value ----------------------------->
# Data table, lists the v and j genes for each chain of each clonotype.
# Function:

get.vdj.info <- function(obj.of.int){
  # ---> Define files.
  cells.clons.info.file <- paste0(gen.data.path, '/FilteredContigAnns_', obj.of.int, '.csv')
  # ---> Load data.
  cells.clons.info <- read.csv(file=cells.clons.info.file, stringsAsFactors=FALSE)
  # ---> Preprocess according to the rules below.
  # Filter out non-productive or not-of-interest contigs from the cells-clonotypes info.
  # Conditions:
  # * For barcodes called as cells.
  # * For clonotype contigs marked with high confidence.
  # * For chains ultimately defined as TRA or TRB.
  # * For productive clonotypes (see 10X support website for a broader definition of 'productive').
  # Also, we'll take out this info since it has been taken into consideration already.
  cells.to.keep <- cells.clons.info$is_cell=='True' & cells.clons.info$high_confidence=='True' & (cells.clons.info$chain=='TRA' | cells.clons.info$chain=='TRB') & cells.clons.info$productive=='True'
  feats.to.keep <- c('barcode', 'contig_id', 'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3', 'cdr3_nt', 'reads', 'umis', 'raw_clonotype_id')
  cells.clons.info <- cells.clons.info[cells.to.keep, feats.to.keep]
  cells.clons.info$clonotype.tag <- cells.clons.info$raw_clonotype_id; cells.clons.info$raw_clonotype_id <- NULL
  cells.clons.info <- as.data.table(cells.clons.info)
  # ---> vdj gene usage info.
  # General vdj gene usage info.
  # We obtain unique v and j genes per clonotype. When cellranger has defined multiple v and j genes for a given clonotype, we take the one that is most supported. The most supported is the gene that has the largest amount of entries (contigs) in the table below.
  vdj.gene.info <- cells.clons.info[
    ,
    .(
      v.gene=.SD[,
        .(count=.N),
        by=v_gene
      ][count==max(count), paste0(unique(v_gene), collapse='|')],
      j.gene=.SD[,
        .(count=.N),
        by=j_gene
      ][count==max(count), paste0(unique(j_gene), collapse='|')]
    ),
    by=.(
      clonotype.tag,
      chain
    )
  ]
  # For those clonotypes with multiple v and j genes that were equally supported, make a guess and keep only one of them.
  vdj.gene.info[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
  vdj.gene.info[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
  # Spread values according to clonotype ID (i.e., to disregard chain information)
  vdj.gene.info[, tmp.genes:=paste(v.gene, j.gene, sep=',')]
  vdj.gene.info <- spread(data=vdj.gene.info[, .(clonotype.tag, chain, tmp.genes)], key='chain', value='tmp.genes', fill=NA)
  # Separate values according to gene type for each chain.
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRA, into=c('tra.v', 'tra.j'), sep=',', convert=FALSE)
  vdj.gene.info <- separate(data=vdj.gene.info, col=TRB, into=c('trb.v', 'trb.j'), sep=',', convert=FALSE)
  # Sanity check. Confirm we have unique entries for the amount of unique clonotypes.
  tmp.check <- vdj.gene.info[, uniqueN(clonotype.tag)==.N]
  if(!tmp.check) stop('Unexpected error post vdj data retrieval.\n')
  # ---> Return processed objects.
  return(vdj.gene.info)
}

# 23 ------------------------------------------------------------------->
# Name: Get table listing the information associated to groups on the basis of clonotype IDs.
# Description:
# Provided the ID of an object, this function provides a table listing the information associated to groups on the basis of clonotype IDs. Clonotypes are defined not only for unique clonotype IDs but also according to donor ID. Then, based on the tag of interest, the table also lists the counts of each group defined for such a tag per clonotype.
# Arguments ------------------------->
# obj.of.int - ID of the object of interest.
# tag.of.int - Tag whose clonotype-specific cell counts' groups will be spread along columns.
# vdj.gene.info - Data table, lists the vdj gene usage information as provided by the function 'get.vdj.info'.
# include.reactivities - Logical, indicates whether or not to include the antigen-reactivity information inferred based on a separate bulk dataset - NOTE: This is quite specific for this project, so you may get rid of all lines related to this option to generalize its use.
# phase.info - Logical, indicats whether the blood draw phase information is defined for the object of interest or if this information should be considered at all.
# cr.info - Logical, indicates whether the cross-reactivity information per clonotype is defined for the object of interest or if this information should be considered at all.
# Value ----------------------------->
# Data table, comprehensive TCR information.
# Function:

get.tcr.by.tag.table <- function(obj.of.int, tmp.tag, vdj.gene.info, include.reactivities=TRUE, phase.info=TRUE, cr.info=TRUE){
  # ---> Prepare data.
  # Retrieve metadata from object of interest.
  meta.data <- cbind(
    data.table(cell=Cells(srt.objs.list[[obj.of.int]])),
    srt.objs.list[[obj.of.int]]@meta.data,
    srt.objs.list[[obj.of.int]]@reductions$umap@cell.embeddings
  )
  # Specify tag of interest (and its respective groups) and add simulated values for tags that may be or not provided.
  if('tmp.tag' %in% colnames(meta.data)) meta.data[, tmp.tag:=NULL]
  meta.data[, tmp.tag:=get(tmp.tag)]
  meta.data[, tmp.tag:=as.character(tmp.tag)]
  tag.groups <- meta.data[, unique(tmp.tag)]
  if(!phase.info | !cr.info){
    if(!phase.info){
      meta.data[, blood.draw.phase.tag:='To be removed']
    }
    if(!cr.info){
      meta.data[, pr.tag:='To be removed']
    }
  }
  # Get info per clonotype.
  # vdj.gene.info <- get.vdj.info(obj.of.int=obj.of.int) # vdj gene usage info. NOTE: It's more practical to have this calculated outside.
  tmp.check <- meta.data[!is.na(clonotype.tag), all(clonotype.tag %in% vdj.gene.info[, clonotype.tag])] # Sanity check. All clonotypes in GEx data should have associated vdj gene usage data.
  if(!tmp.check) stop('Unexpected error 1.\n')
  clones.info <- meta.data[
    !is.na(clonotype.tag),
    .(
      trb.nt=unique(TRB.nt.chains.tag),
      tra.nt=unique(TRA.nt.chains.tag),
      trb.aa=unique(TRB.aa.chains.tag),
      tra.aa=unique(TRA.aa.chains.tag),
      global.size=unique(clon.size.tag)
    ),
    by=.(clonotype=clonotype.tag)
  ]
  clones.info <- merge(x=clones.info, y=vdj.gene.info, by.x='clonotype', by.y='clonotype.tag', all.x=TRUE, all.y=FALSE)
  # Get cell counts per donor, per phase, per value of tag of interest.
  # NOTE: You may get rid of the first chunk within (and thus just keep the else section) if you wish to generalize the function use.
  if(include.reactivities){
    tmp.data <- meta.data[
      !is.na(clonotype.tag) & !is.na(donor.id.tag),
      .(
        cell.count=.N,
        blood.flu=`TCR-L.FLU`,
        blood.sars=`TCR-L.SARS`,
        blood.cr=`TCR-L.SARS-FLU`,
        consensus.reactivity=reactivity.tag
      ),
      by=.(
        clonotype=clonotype.tag, pr.tag,
        donor=donor.id.tag, phase=blood.draw.phase.tag,
        tmp.tag
      )
    ]
    tmp.data[blood.flu=='Not', blood.flu:=NA]
    tmp.data[blood.sars=='Not', blood.sars:=NA]
    tmp.data[blood.cr=='Not', blood.cr:=NA]
  }else{
    tmp.data <- meta.data[
      !is.na(clonotype.tag) & !is.na(donor.id.tag),
      .(
        cell.count=.N,
        consensus.reactivity=reactivity.tag
      ),
      by=.(
        clonotype=clonotype.tag, pr.tag,
        donor=donor.id.tag, phase=blood.draw.phase.tag,
        tmp.tag
      )
    ]
  }
  tmp.data <- unique(tmp.data)
  # Spread data according to the groups defined for the tag of interest.
  tmp.data <- spread(data=tmp.data, key=tmp.tag, value=cell.count, fill=0)
  # Include local size.
  tmp.vals <- rowSums(as.matrix(tmp.data[, ..tag.groups]))
  tmp.data[, local.size:=tmp.vals]
  # Merge counts with clonotype data.
  tmp.data <- merge(x=clones.info, y=tmp.data, by='clonotype')
  # Formatting.
  cols.order <- c('clonotype', 'trb.aa', 'tra.aa', 'trb.nt', 'tra.nt', 'trb.v', 'tra.v', 'trb.j', 'tra.j', 'pr.tag', 'donor', 'phase', 'global.size', 'local.size', mixedsort(tag.groups))
  if(include.reactivities) cols.order <- c(cols.order, c('blood.sars', 'blood.flu', 'blood.cr', 'consensus.reactivity'))
  tmp.data <- tmp.data[, ..cols.order]
  tmp.data[, to.sort:=as.integer(str_replace(string=clonotype, pattern='clonotype', replacement=''))]
  setorderv(x=tmp.data, cols='to.sort'); tmp.data[, to.sort:=NULL]
  return(tmp.data)
}


# 24 ------------------------------------------------------------------->
# Name: Get cell counts per donor per cluster.
# Description:
# TBD.
# Arguments ------------------------->
# table.obj - ID of the object of interest.
# Value ----------------------------->
# Data table, cell counts per cluster.
# Function:

get.sep.conv.table <- function(table.obj, set.lab){
  # @ General definitions.
  # Current object.
  seurat.obj <- srt.objs.list[[table.obj]]
  clusters.tag <- clusts.lab[table.obj]
  meta.data <- cbind(
    data.table(cell=Cells(seurat.obj)),
    seurat.obj@meta.data,
    seurat.obj@reductions$umap@cell.embeddings
  )
  meta.data[, cluster:=get(clusters.tag)]
  # Previous object.
  #   Define object.
  prev.seurat.obj <- get.raw.obj(obj.of.int=table.obj)
  #   Replace HTO negatives w/ Demuxlet
  prev.seurat.obj <- rep.hto.negatives.with.demuxlet(seurat.obj=prev.seurat.obj, set.lab=set.lab, tissue.type=ifelse(test=str_detect(string=set.lab, pattern='tumor'), yes='Tumor-lung', no='Normal-lung'), aggr.pffx.lab='R24_Cancer_Batches-1-to-20_')
  prev.seurat.obj@meta.data[, 'donor.id.tag'] <- as.integer(prev.seurat.obj@meta.data[, 'full.donor.id.tag'])
  prev.meta.data <- cbind(
    data.table(cell=Cells(prev.seurat.obj)),
    prev.seurat.obj@meta.data
  )
  #   Exclude barcodes classified as doublets based on Hashtag data.
  prev.meta.data <- prev.meta.data[hashtag.tag!='Doublet' | is.na(hashtag.tag)]
  # Get cell counts per groups per cluster for the current object.
  tmp.data <- meta.data[,
    .(cell.count=.N),
    by=.(
      donor=donor.id.tag,
      cluster=cluster
    )
  ]
  clusters <- meta.data[, unique(as.character(cluster))]
  tmp.data <- spread(data=tmp.data, key=cluster, value=cell.count, fill=0)
  # Get cell counts per groups for the current object.
  post.counts <- meta.data[,
    .(post.counts=.N),
    by=.(
      donor=donor.id.tag
    )
  ]
  tmp.data <- merge(x=tmp.data, y=post.counts, by=c('donor'))
  # Get cell counts per groups for the previous object.
  pre.counts <- prev.meta.data[,
    .(pre.counts=.N),
    by=.(
      donor=donor.id.tag
    )
  ]
  tmp.data <- merge(x=tmp.data, y=pre.counts, by=c('donor'))
  tmp.data$donor <- as.character(tmp.data$donor)
  tmp.data[is.na(donor), donor:='Unassigned']
  # @ Sanity checks
  # No entry should have unassigned values except for the indicated datasets (virus tag was assigned as NA).
  to.check <- any(apply(X=as.matrix(tmp.data), MARGIN=1, FUN=function(x) any(is.na(x))))
  if(to.check) stop('Unexpected error 1.\n')
  # Pre-filtering counts should be larger than the post-filtering counts.
  # to.check <- tmp.data[, all(pre.counts>=post.counts)]
  # if(to.check) stop('Unexpected error 2.\n')
  # @ Final formatting.
  cols.order <- c('donor', 'pre.counts', 'post.counts', clusters)
  names(cols.order) <- c('Donor', 'Pre', 'Post', clusters)
  tmp.data <- tmp.data[, ..cols.order]
  colnames(tmp.data) <- names(cols.order)
  return(tmp.data)
}


# 25 ------------------------------------------------------------------->
# Name: Output UpSet plot from metadata according to tag of interest.
# Description: Provided the metadata from a seurat object and a tag of interest, this script will output an UpSet plot showing the clonotype sharing by the groups under the tag of interest.
# Arguments ------------------------->
#     @ About data.
#         Data specificities.
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# size.thold - Clone size threshold to consider when picking clonotypes to depict in plot.
#     @ Data filtering and grouping.
# groups.tag - Tag defined in the seurat object metadata whose groups will be represent the sets in the UpSet plot. No default.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# filter.tag - Tag defined in the seurat object metadata whose groups will be used to filter out cells. No cell filtering is applied when this argument's set to NULL. Default: NULL.
# groups.to.filter - Character listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
#     @ About plot and output.
# TBD.
#         Output.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

output.upset.plot <- function(
  # @   About data.
  # Data specificities.
  seurat.obj, size.thold=1,
  # Data filtering and grouping.
  groups.tag, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE,
  # @   About plot and output.
  # Color scale
  # this.color.scale=NULL, col.min=NULL, col.max=NULL,
  # Size scale.
  # scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  # File.
  file.name=NULL
){
  # ---> Fetch data.
  # Get clone size.
  feature.vals <- try(expr=as.data.table(FetchData(object=seurat.obj, vars=c('clonotype.tag', 'clon.size.tag'))), silent=TRUE)
  # Kill if it was not possible to find the feature in the data slot appointed.
  if(any(class(feature.vals)=='try-error')) stop('Seurat object does not seem to be annotated. Metadata should include columns "clonotype.tag" and "clon.size.tag". Please provide a validly annotated object.\n')
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)) if(!(filter.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for filtering set to a value different than NULL, but not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)){
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag], filter.tag=seurat.obj@meta.data[, filter.tag])
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='rest']
  }
  # ---> Filter out non-desired cells.
  # According to filter tag.
  if(!is.null(filter.tag) & !is.null(groups.to.filter)){
    if(keep) tmp.data <- tmp.data[filter.tag %in% groups.to.filter] else tmp.data <- tmp.data[!filter.tag %in% groups.to.filter]
    tmp.data[, filter.tag:=NULL]
  }
  # ---> Collect data to be used as UpSet input.
  these.groups <- tmp.data[, unique(groups.tag)]
  upset.input <- lapply(X=these.groups, FUN=function(this.group){
    tmp.data <- tmp.data[groups.tag==this.group & !is.na(clonotype.tag), .(cell.count=.N), by=.(clone=clonotype.tag)]
    tmp.data <- tmp.data[cell.count>=size.thold, clone]
    return(tmp.data)
  })
  names(upset.input) <- these.groups
  # Output.
  # Normal version.
  tmp.ggplot <- upset(
    data=fromList(upset.input), main.bar.color='#0099cc'
  )
  tmp.file.name <- paste0(file.name, '.C.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
  # Blank version
  tmp.ggplot <- upset(
    data=fromList(upset.input), main.bar.color='#0099cc', text.scale=0
  )
  tmp.file.name <- paste0(file.name, '.B.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
  # Mock.
  return(NA)
}


# 26 ------------------------------------------------------------------->
# Name: List samples' information.
# Description: The following function retrieves the information provided in the original aggr and vdj aggr tables from a set of aggregations of interest and then returns such information.
# TBD.
# Arguments ------------------------->
# table.obj - ID of the object of interest.
# Value ----------------------------->
# Data table, aggr tables' information.
# Function:

list.samples.info <- function(table.obj){
  # Load aggr table.
  aggr.table.file <- paste0(gen.data.path, '/ForAnns_', table.obj, '.csv')
  if(!file.exists(aggr.table.file)) stop(paste0('Aggregation table file does not exist.\n', aggr.table.file, '\n'))
  vdj.table.file <- paste0(gen.data.path, '/ForAnns_vdj_', table.obj, '.csv')
  if(!file.exists(vdj.table.file)) warning(paste0('Aggregation vdj table file does not exist.\n', vdj.table.file, '\n'))
  aggr.table <- fread(file=aggr.table.file)
  # Load vdj aggr table only if necessary.
  if(file.exists(vdj.table.file)){
    vdj.table <- fread(file=vdj.table.file)
    # Sanity checks.
    tmp.check.1 <- nrow(aggr.table)==nrow(vdj.table)
    tmp.check.2 <- all(aggr.table[, library_id]==vdj.table[, library_id])
    tmp.check <- tmp.check.1 & tmp.check.2
    if(!tmp.check) stop('Unexpected error.\n')
  }else{
    vdj.table <- NULL
  }
  # Get unique IDs represented by the whole path to molecule h5 file.
  aggr.table[, sample.id:=molecule_h5]
  aggr.table[, sample.id:=str_replace(string=sample.id, pattern='/molecule_info.h5$', replacement='')]
  # Retrieve vdj path
  if(!is.null(vdj.table)) aggr.table[, vdj.path:=str_replace(string=vdj.table[, clonotypes], pattern='/clonotypes.csv', replacement='')]
  # Retrieve lane tag information.
  lane.info <- separate(data=aggr.table[, .(lane.tag)], col=lane.tag, into=str_split(string=lane.tags[table.obj], pattern=';', simplify=TRUE)[1, ], sep=';')
  aggr.table <- cbind(aggr.table, lane.info)
  return(aggr.table)
}


# 27 ------------------------------------------------------------------->
# Name: Get comprehensive reactivity summary
# Description:
# Obtain raw numbers and fractions of ag-reactive items per donor, where item may be cell or clonotype.
# Arguments ------------------------->
# obj.of.int - ID of the object of interest.
# Value ----------------------------->
# Data table, raw numbers and fractions of ag-reactive items per donor.
# Function:

get.comp.reactivity.summ <- function(obj.of.int, cell.count=TRUE){
  # @ General definitions.
  # Current object.
  seurat.obj <- srt.objs.list[[obj.of.int]]
  clusters.tag <- clusts.lab[obj.of.int]
  meta.data <- cbind(
    data.table(cell=Cells(seurat.obj)),
    seurat.obj@meta.data,
    seurat.obj@reductions$umap@cell.embeddings
  )
  meta.data[, cluster:=get(clusters.tag)]
  meta.data[is.na(reactivity.tag) & !is.na(clonotype.tag), reactivity.tag:='Unknown']
  # ---> Obtain general data.
  if(cell.count){
    # For cell-specific summary.
    tmp.data.1 <- meta.data[
      !is.na(donor.id.tag) & !is.na(clonotype.tag),
      .(ind.count=.N),
      by=.(donor=donor.id.tag, reactivity=reactivity.tag)
    ]
    tmp.data.2 <- meta.data[
      !is.na(donor.id.tag) & !is.na(clonotype.tag),
      .(total.count=.N),
      by=.(donor=donor.id.tag)
    ]
  }else{
    # For clonotype-specific summary.
    tmp.data.1 <- meta.data[
      !is.na(donor.id.tag) & !is.na(clonotype.tag),
      .(ind.count=uniqueN(clonotype.tag)),
      by=.(donor=donor.id.tag, reactivity=reactivity.tag)
    ]
    tmp.data.2 <- meta.data[
      !is.na(donor.id.tag) & !is.na(clonotype.tag),
      .(total.count=uniqueN(clonotype.tag)),
      by=.(donor=donor.id.tag)
    ]
  }
  # ---> Obtain raw numbers (tidy table).
  to.output.1 <- spread(data=tmp.data.1, key=reactivity, value=ind.count, fill=0)
  # ---> Obtain proportions (tidy table).
  tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor', all=TRUE)
  tmp.data[, ind.prop:=ind.count/total.count]
  tmp.data <- tmp.data[, .(donor, reactivity, ind.prop)]
  to.output.2 <- spread(data=tmp.data, key=reactivity, value=ind.prop, fill=0)
  # ---> Merge results and output.
  to.output <- merge(x=to.output.1, y=to.output.2, by='donor', suffixes=c('.C', '.P'))
  return(to.output)
}


# 28 ------------------------------------------------------------------->
# Name: Harmonize FACS-sorting batch (for the DICE Cancer project)
# NOTE: THIS IS QUITE A SPECIFIC STEP FOR THE DICE CANCER PROJECT
# NOTE: THIS IS QUITE A SPECIFIC STEP FOR THE HEALTHY TISSUE SUBPROJECT.
# Description:
# Harmonize FACS sorting batch. The batches for the healthy adjacent tissue are labeled by tetters that are indexed by an even number, whilst batches for tumor tissue are labeled by letters that are indexed by odd numbers. Thus, we must turn them into a sequential list of letters. For example, B -> A and D -> C. Of NOTE, this does not apply to FACS-sorting batch W.
# Arguments ------------------------->
# batch.vector - Charecter vector, represents batches letters.
# Value ----------------------------->
# Harmonized vector of batch IDs.
# Function:

harmonize.facs.batches <- function(batch.vector){
  # Create reference
  batch.ref <- LETTERS[1:length(unique(batch.vector))]
  names(batch.ref) <- sort(unique(batch.vector))
  # Translate accordingly.
  batch.vector <- batch.ref[batch.vector]
  batch.vector <- unname(batch.vector)
  return(batch.vector)
}


# 29 ------------------------------------------------------------------->
# Name: Harmonize sequencing batch (for the DICE Cancer project)
# NOTE: THIS IS QUITE A SPECIFIC STEP FOR THE DICE CANCER PROJECT
# NOTE: THIS IS QUITE A SPECIFIC STEP FOR THE HEALTHY TISSUE SUBPROJECT.
# Description:
# Harmonize sequencing batch. Set a specific order for the sequencing batches out the general sequencing batch that was assigned when taking into consideration also the tumor tissue libraries.
# Arguments ------------------------->
# batch.vector - Charecter vector, represents sequencing batches.
# Value ----------------------------->
# Harmonized vector of batch IDs.
# Function:

harmonize.seq.batches <- function(batch.vector){
  # Create reference
  uniq.entries <- mixedsort(unique(unlist(str_split(string=batch.vector, pattern=';', simplify=FALSE))))
  batch.ref <- 1:length(uniq.entries)
  names(batch.ref) <- uniq.entries
  # Translate accordingly.
  batch.vector <- lapply(X=batch.vector, FUN=function(tmp.batches){
    tmp.batches <- str_split(string=tmp.batches, pattern=';')[[1]]
    if(!all(tmp.batches %in% names(batch.ref))) stop('Unexpected error 1 in function "harmonize.seq.batches".\n')
    tmp.batches <- paste0(batch.ref[tmp.batches], collapse=';')
    return(tmp.batches)
  })
  batch.vector <- unlist(batch.vector)
  return(batch.vector)
}


# 30 ------------------------------------------------------------------->
# Name: Harmonize library ID (for the DICE Cancer project)
# NOTE: THIS IS QUITE A SPECIFIC STEP FOR THE DICA CANCER PROJECT
# Description:
# Harmonize library ID as it was done originally when defining the library IDs in the aggr tables.
# Arguments ------------------------->
# sample.sheet - Datatable, represents a sample sheet where the column "library_id" must be defined.
# Value ----------------------------->
# Datatable with the same format as the input except for the column "library_id" which is harmonized here.
# Function:

harmonize.ss <- function(sample.sheet){
  # Donors' sorting batch, first part.
  sorting.batch.1 <- str_extract(string=sample.sheet$library_id, pattern='_Hu_\\d+[^_]*')
  sorting.batch.1 <- str_replace(string=str_replace(string=sorting.batch.1, pattern='_Hu_', replacement=''), pattern='\\D+', replacement='')
  # Donors' sorting batch, second part.
  # @ For the first pattern available.
  sorting.batch.1.tmp <- str_extract(string=sample.sheet$library_id, pattern='_[A-Z]{1}_\\d+D$')
  sorting.batch.1.tmp <- str_replace(string=str_replace(string=sorting.batch.1.tmp, pattern='_\\d+D$', replacement=''), pattern='_', replacement='')
  # @ For the second pattern available.
  sorting.batch.2.tmp <- str_extract(string=sample.sheet$library_id, pattern='_Hu_\\d+[^_]*')
  sorting.batch.2.tmp <- str_replace(string=sorting.batch.2.tmp, pattern='_Hu_\\d+', replacement='')
  sorting.batch.2.tmp[sorting.batch.2.tmp==''] <- NA
  # @ Coalesce batches from all patterns and replace NA values with 'A' (i.e., indicating there was a single batch).
  sorting.batch.2 <- dplyr::coalesce(sorting.batch.1.tmp, sorting.batch.2.tmp)
  sorting.batch.2[is.na(sorting.batch.2)] <- 'A'
  # Combine both parts of the donors' sorting batch.
  sorting.batch <- paste0(sorting.batch.1, sorting.batch.2)
  # Get number of donors.
  number.of.donors <- str_extract(string=sample.sheet$library_id, pattern='\\d+D$')
  # Add appropriate sample ID by eliminating redundant information.
  sample.sheet$library_id <- str_replace(string=sample.sheet$library_id, pattern='(_[A-Z]{1})*_\\d+D$', replacement=paste0('_', sorting.batch, '_', number.of.donors))
  sample.sheet$library_id <- str_replace(string=sample.sheet$library_id, pattern='_Hu_[^_]+', replacement='_Hu')
  # Extract FACS-sorting batch info.
  tmp.batch.info <- str_extract(string=sample.sheet$library_id, pattern='^R24[^_]+|R24_Cancer')
  tmp.batch.info <- str_replace(string=tmp.batch.info, pattern='R24', replacement='')
  tmp.batch.info[tmp.batch.info=='_Cancer'] <- 'A'
  # Modify library name according to the name of batch A.
  sample.sheet$library_id <- str_replace(string=sample.sheet$library_id, pattern='^R24[^_]+|R24_Cancer', replacement=paste0('R24-Cancer-', tmp.batch.info))
  # Return
  return(sample.sheet)
}
