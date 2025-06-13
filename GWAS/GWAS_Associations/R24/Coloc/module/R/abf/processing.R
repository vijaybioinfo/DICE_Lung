library(GenomicRanges)
source("utils.R")

##############################################################################################################################################
############################################################ DATA PROCESSING FUNCTIONS #######################################################
##############################################################################################################################################

#######################################################
################# CASE IDENTIFICATION #################
#######################################################
identify.case <- function(data){
    cols <- colnames(data)

    if ((COLS$BETA %in% cols) & (COLS$SE %in% cols)) return(1)
	if ((COLS$PVAL %in% cols) & (COLS$AF %in% cols)) return(2)

    stop("Case could not be identified. You must provide either (BETA and SE) or (P and AF)")
    
}


#######################################################
################ CHECK DATA REQUIREMENTS ##############
#######################################################
check.case.1.requirements <- function(data, type, cc=NULL, N=NULL){
	cols <- colnames(data)

    if (!(COLS$BETA %in% cols) || !(COLS$SE %in% cols)) stop("Case 1 requires BETA and SE columns")
		
    if (type == "quant"){
        if (!(COLS$AF %in% cols) || is.null(N)) stop("Case 1 and type 'quant' require allele frequency (AF) and sample size")
    }
}

check.case.2.requirements <- function(data, type, cc=NULL, N=NULL){
	cols <- colnames(data)

    if (!(COLS$PVAL %in% cols) || !(COLS$AF %in% cols) || is.null(N)) stop("Case 2 requires pvalue(P) and allele frequency(AF) columns, and sample size")

	if (type == "cc" & is.null(cc)) stop("Case 2 and type 'cc' requires case control ratio")
}



#######################################################
################### HANDLE CASES ######################
#######################################################
prepare.for.case.1 <- function(data, type="quant", cc=0.3, N=NULL){
    ### Check case requirements ###
	check.case.1.requirements(data, type, cc, N)

	### Remove problematic values ###
	if (type == "quant"){
		core.columns <- list(CHR = COLS$CHR, POS = COLS$POS, BETA =  COLS$BETA, SE = COLS$SE, AF = COLS$AF)
	} else {
		core.columns <- list(CHR = COLS$CHR, POS = COLS$POS, BETA = COLS$BETA, SE = COLS$SE)
	}
	
	data <- data[rowSums(is.na(data[, unlist(core.columns)])) == 0, ]
	data <- data[rowSums(data[, unlist(core.columns)] == Inf) == 0, ]

	### Estimate VARBETA ###
	data$varbeta <- data[[COLS$SE]] ^ 2
	data <- data[data$varbeta != 0, ]

	if (type == "quant"){
		data <- data[(data[[COLS$AF]] > 0) & (data[[COLS$AF]] < 1), ]
		
		return (list(data=data, meta=list(case=1, type=type, N=N)))
	}
	
	return (list(data=data, meta=list(case=1, type=type, s=cc)))
}



prepare.for.case.2 <- function(data, type="quant", cc=0.3, N=NULL){
    ### Check case requirements ###
	check.case.2.requirements(data, type, cc, N)

    ### Remove problematic values ###
	core.columns <- list(CHR=COLS$CHR, POS=COLS$POS, P=COLS$PVAL, AF=COLS$AF)

	data <- data[rowSums(is.na(data[, unlist(core.columns)])) == 0, ]
    data <- data[rowSums(data[, unlist(core.columns)] == Inf) == 0, ]

    data <- data[data[[core.columns$P]] > 0, ]
    data <- data[(data[[core.columns$AF]] > 0) & (data[[core.columns$AF]] < 1), ]
	
    if (type == "cc"){
        return (list(data=data, meta=list(case=2, type=type, s=cc, N=N)))
    }

    return (list(data=data, meta=list(case=2, type=type, N=N)))
}


CASE.HANDLERS <- list(prepare.for.case.1, prepare.for.case.2)



########################################################
################# PREPARE DATA OBJECT ##################
########################################################
prepare.data.objects <- function(data, type="quant", cc=0.3, N=NULL){
	### Identify case ###
	case <- identify.case(data)

    ### Prepare ###
    prepare.case <- CASE.HANDLERS[[case]]
    
    return (prepare.case(data, type, cc, N))
}





#######################################################################################################################################################
############################################################# MAKE DATASETS ###########################################################################
#######################################################################################################################################################

#######################################################
############# MERGE SNPS FUNCTION BUNDLE ##############
#######################################################
.to.ranges <- function(data){
	return(GRanges(seqnames=data[[COLS$CHR]], ranges=IRanges(start=data[[COLS$POS]], width=1)))
}

.overlap <- function(d1, d2){
	print("Overlapping...")
	gr1 <- .to.ranges(d1)
	gr2 <- .to.ranges(d2)
	
	overlap <- findOverlaps(gr1, gr2)
	d1 <- d1[queryHits(overlap), ]
	d2 <- d2[subjectHits(overlap), ]

	return( list(dataset1=d1, dataset2=d2) )
}

.make.pairs <- function(d1, d2){
	print("Making pairs...")

	pairs <- paste0(d1[,"Feature1"], "_", d2[,"Feature2"])

	d1["Pair"] <- pairs
	d2["Pair"] <- pairs
	
	return( list(dataset1=split(d1, d1$Pair), dataset2=split(d2, d2$Pair)) )
}

.keep.intersecting.snps <- function(d1, d2){
	print("Merging SNPS...")

	d1 <- d1[!duplicated(d1[[COLS$SNPID]]), ]
	d2 <- d2[!duplicated(d2[[COLS$SNPID]]), ]
	
	snpid_intersection <- intersect(d1[[COLS$SNPID]], d2[[COLS$SNPID]])
	d1 <- d1[d1[[COLS$SNPID]] %in% snpid_intersection, ]
	d2 <- d2[d2[[COLS$SNPID]] %in% snpid_intersection, ]
	
	d1 <- d1[order(d1[[COLS$POS]]), ]
	d2 <- d2[order(d1[[COLS$POS]]), ]
	
	return(list(d1=d1, d2=d2))
}

.merge.snps <- function(d1, d2){
	########### Overlap datasets ###############
	shared <- .overlap(d1, d2)
	
	if (check.empty(shared$dataset1) | check.empty(shared$dataset2)){
		print("At least 1 of the dataset is empty...")
		return(NULL)
	}

	########### Make pairs of entities #############
	pairs <- .make.pairs(shared$dataset1, shared$dataset2)
	d1 <- pairs$dataset1
	d2 <- pairs$dataset2

	########### Keep only snps that are shared ##############
	intersecting <- mapply(.keep.intersecting.snps, d1=d1, d2=d2, SIMPLIFY = FALSE)

	indexes <- which(sapply(intersecting, function(X){!(check.empty(X$d1) | check.empty(X$d2))}))
	intersecting <- intersecting[names(intersecting)[indexes]]
	
	if (length(intersecting) == 0){
		return(NULL)
	}

	return (intersecting)
}



########################################################
################# MAKE COLOC DATASETS ##################
########################################################
make.datasets <- function(D1, D2){
	########## Merge SNPS ###########
	objects <- .merge.snps(D1$data, D2$data)

	if (is.null(objects)) return(NULL)

	########## Add metadata #########
	datasets <- lapply(objects, function(O){
							list(
								d1 = c(as.list(keep.req.columns(O$d1, D1$meta$case, D1$meta$type)), D1$meta),
								d2 = c(as.list(keep.req.columns(O$d2, D2$meta$case, D2$meta$type)), D2$meta),
								tables = list(d1 = O$d1, d2 = O$d2)
							)
						}
					)
	return(datasets)	
}