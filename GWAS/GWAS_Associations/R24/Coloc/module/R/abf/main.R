# libraries specifically for colocalization
library(coloc)
library(purrr)

source('utils.R')
source('summary.R')
source("processing.R")

options(scipen = 10)


########################
######## TODO ##########
########################

# 1.- Source each file into a different env instead of loading everything into the glob space

########### Other Functions ##############
kill <- function(args){
	empty.objects <- make.empty()
    write.coloc.output(empty.objects$top, empty.objects$full, args$out)
	print("Killing because empty objects")
    quit(status=0)
}

kill.if.empty <- function(data, args){
    if (check.empty(data)) kill(args)
}

kill.if.null <- function(data, args){
	if (is.null(data)) kill(args)
}

##########################################################################################
################################ COLOC ABF METHOD ########################################
##########################################################################################
run.abf <- function(datasets, sensitivity = NULL){
	###### Run Colocalization ########
	results <- lapply(datasets, FUN=function(X){ coloc.abf(dataset1=X$d1, dataset2=X$d2) })

	###### Create Summary Tables #####
	summarized <- lapply(names(results), FUN=function(X){ summarize(results[[X]], datasets[[X]])})
	top.associations <- rbindlist(lapply(summarized, FUN=function(X){X$top}), fill=TRUE)
	full.associations <- rbindlist(lapply(summarized, FUN=function(X){X$full}), fill=TRUE)
	
	###### Make Sensitivity plots ######
	if (!(sensitivity == FALSE | is.null(sensitivity))){
		sensitivity <- paste0(sensitivity, "/sensitivity/")
        create.ouput.folder(sensitivity)
		lapply(names(results), FUN=function(name){ plot.sensitivity(results[[name]], paste0(sensitivity, name, ".pdf")) })
	}

	return(list(top=top.associations, full=full.associations))
}

#################################################################################################
################################# MAIN WRAPPER FOR run.abf ######################################
#################################################################################################
main.abf <- function(args){
	########### read datasets ################
	d1 <- read(args$d1, args$feature1, COLS$F1)
	d2 <- read(args$d2, args$feature2, COLS$F2)

	########### Prepare data objects ###########
	obj1 <- prepare.data.objects(d1, args$type1, args$cc_ratio1, args$N1)
	obj2 <- prepare.data.objects(d2, args$type2, args$cc_ratio2, args$N2)

	kill.if.empty(obj1$data, args)
	kill.if.empty(obj2$data, args)

	########### Make coloc datasets #############
	datasets <- make.datasets(obj1, obj2)
	kill.if.null(datasets, args)

	########### Remove objects ###############
	rm(d1)
	rm(d2)

	########### Run analysis #############
	results = run.abf(datasets, sensitivity=args$out)

	return(results)
}

###########################################################################################
######################### RUN COLOC ABF FROM COMMAND LINE #################################
###########################################################################################
cli <- function(arguments){
    library(argparse)

    parser <- ArgumentParser()

	parser$add_argument("--d1", type="character", help="Summary statistics file 1", required=TRUE)
	parser$add_argument("--d2", type="character", help="Summary statistics file 2", required=TRUE)
	parser$add_argument("--feature1", default="Loci", type="character", help="Column name for molecular feature on dataset 1 (e.g. GWAS locus, eGene, chromatin peak)")
	parser$add_argument("--feature2", default="GENEID", type="character", help="Column name for molecular feature on dataset 2 (e.g. GWAS locus, eGene, chromatin peak)")
	parser$add_argument("--type1", default="quant", help="Type of trait 1 {quant, cc}")
	parser$add_argument("--type2", default="quant", help="Type of trait 2 {quant, cc}")
	parser$add_argument("--N1", type="numeric", default=NULL, help="Sample size of trait 1 ")
	parser$add_argument("--N2", type="numeric", default=NULL, help="Sample size of trait 2 ")
	parser$add_argument("--cc-ratio1", type="double", default=0.3, help="If type1 is 'cc' provide case-control ratio")
	parser$add_argument("--cc-ratio2", type="double", default=0.3, help="If type2 is 'cc' provide case-control ratio")
	parser$add_argument("--out", type="character", help="Output file prefix", required=TRUE)
	
    args <- parser$parse_args(arguments)
    return(args)
}

args.default.list <- list(
						d1 = NULL,
						d2 = NULL,
						feature1 = "Loci",
						feature2 = "Loci",
						type1 = "quant",
						type2 = "quant",
						N1 = NULL,
						N2 = NULL,
						cc_ratio1 = 0.3,
						cc_ratio2 = 0.3,
						out = NULL
						)

run.from.cli <- function(arguments){
	########### Load CLI ##############
	args <- cli(arguments)

	########### Run Main analysis ###############
	results <- main.abf(args)

	################# Write results #######################
	write.coloc.output(results$top, results$full, R.utils::getAbsolutePath(args$out))
}

if (sys.nframe() == 0L){
	run.from.cli(commandArgs(trailingOnly = TRUE))
}




