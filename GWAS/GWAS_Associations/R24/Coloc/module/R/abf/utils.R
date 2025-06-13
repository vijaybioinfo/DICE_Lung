library(data.table)

options(datatable.fread.datatable=FALSE)

###########################################################################################################################
###################################################### COLUMNS ############################################################
###########################################################################################################################

COLS <- list(
			CHR="Chromosome", 
			POS="Position", 
			SNPID="SNPID", 
			PVAL="P", 
			BETA="BETA", 
			SE="SE", 
			AF="AF", 
			N="N", 
			VARBETA="varbeta",
			F1 = "Feature1",
			F2 = "Feature2",
			PPH0 = "PPH0",
			PPH1 = "PPH1",
			PPH2 = "PPH2",
			PPH3 = "PPH3",
			PPH4 = "PPH4",
			PPH4.PPH3 = "PPH4/PPH3",
			SNP.PP.H4 = "SNP.PP.H4",
			METHOD = "METHOD"
			)



########### METHODS ##########
keep.req.columns <- function(data, case, type){
	cols <- switch(case,
					switch(type, "quant" = c(COLS$SNPID, COLS$BETA, COLS$VARBETA, COLS$AF), "cc" = c(COLS$SNPID, COLS$BETA, COLS$VARBETA)),
					c(COLS$SNPID, COLS$PVAL, COLS$AF)
				)

	data <- data[, cols]

	TRANSLATOR <- list("SNPID"="snp", "BETA"="beta", "varbeta"="varbeta", "P"="pvalues", "AF"="MAF")
	indexes <- match(names(TRANSLATOR), names(data))

	names(data)[indexes[!is.na(indexes)]] <- unlist(TRANSLATOR[!is.na(indexes)])
	
	return(data)
}


##########################################################################################################################
####################################################### HELPER FUNCTIONS #################################################
##########################################################################################################################
check.empty <- function(data, args){
	return(nrow(data) == 0)
}

make.empty <- function(){
	############## Define Columns ################
	core.columns <- c(COLS$CHR, COLS$F1, COLS$F2, COLS$PPH0, COLS$PPH1, COLS$PPH2, COLS$PPH3, COLS$PPH4, COLS$PPH4.PPH3, COLS$SNPID, COLS$SNP.PP.H4, COLS$POS)
	stats.columns <- c(COLS$PVAL, COLS$BETA, COLS$SE, COLS$AF)

	top <- data.frame(matrix(nrow=0, ncol=13))
	colnames(top) <- c(core.columns, COLS$METHOD)
	full <- data.frame(matrix(nrow=0, ncol=21))
	colnames(full) <- c(core.columns, paste0("F1_", stats.columns), paste0("F2_", stats.columns), COLS$METHOD)

	return(list(top=top, full=full))
}

###############################################################################################################
################################################# IO FUNCTIONS ################################################
###############################################################################################################

###########################
########## INPUT ##########
###########################

read <- function(file, feature=NULL, rename=NULL){
    print("Reading...")
    data <- fread(file)
    data[data[[COLS$CHR]] == "True", COLS$CHR] <- 1

	if ("MAF" %in% colnames(data)){
		data[[COLS$AF]] <- data$MAF
	}

    ### Final touch ###
    if(!is.null(feature)){
		if (!(feature %in% colnames(data))){
			stop(paste0("Feature ", feature, " not in columns"))
		}

		if (is.null(rename)) stop("New column name for feature not provided")
		colnames(data)[colnames(data) == feature] <- rename
	}

    return(data)
}

###########################
########## OUTPUT #########
###########################

write.output <- function(object, filename){
	fwrite(object, filename, sep="\t", row.names=FALSE, quote=FALSE)
}

create.ouput.folder <- function(path){
	if (!file.exists(path)) dir.create(path, recursive=TRUE)
}

write.coloc.output <- function(association, all, prefix){
    create.ouput.folder(dirname(prefix))
    write.output(association, paste0(prefix,".tsv"))
    write.output(all, paste0(prefix,"_snps.tsv"))
}


###########################################################################################################################################
############################################################## SENSITIVITY ANALYSIS #######################################################
###########################################################################################################################################
plot.sensitivity <- function(results, out){
	########## Make sensitivity plot #########
	if(!(is.nan(results$summary["PP.H4.abf"])) & results$summary["PP.H4.abf"] >= 0.5){

		pdf(out)
		try(sensitivity(results, rule="H4 > 0.5"))
		dev.off()

	} else {
		print("Sensitivity plot can not be created... skiping...")
	}
}



