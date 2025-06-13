source("utils.R")

### TODO: FINISH CODE
.add.prefix <- function(x, prefix){
	return(paste0(prefix, "_", x))
}

.set_columns <- function(data, prefix){
	Feature <- grep("Feature", names(data), value = TRUE)

	if (length(Feature) != 1){
		print(names(data))
		stop("There must be one single Feature column...")
	}

	columns <- c(COLS$CHR, COLS$POS, COLS$SNPID, Feature, COLS$P, COLS$BETA, COLS$SE, COLS$AF)
	no.exists <- columns[!(columns %in% names(data))]

	d <- data
	d[, no.exists] <- NA
	d <- d[, columns]
	colnames(d) <-  c(COLS$CHR, COLS$POS, COLS$SNPID, Feature, paste0(prefix, "_", c(COLS$P, COLS$BETA, COLS$SE, COLS$AF)))
	
	return(d)
}

summarize <- function(results, datasets){
	if (is.null(results)){
		return(make.empty())
	}
	
	############ SNP-level summary #############
	snps <- as.data.frame(results$results[, c("snp", "SNP.PP.H4")])
	colnames(snps) <- c(COLS$SNPID, COLS$SNP.PP.H4)
	
	############ Retrieve summary stats columns ####################
	d1 <- .set_columns(datasets$tables$d1, "F1")
	d2 <- .set_columns(datasets$tables$d2, "F2")
	d2[, c(COLS$CHR, COLS$POS)] <- NULL

	############ Merge INFO ################
	merged <- merge(d1, snps, by=COLS$SNPID)
	merged <- merge(merged, d2, by=COLS$SNPID)

	############ Get Posterior Probabilities ################
	merged[[COLS$PPH0]] <- results$summary["PP.H0.abf"]
	merged[[COLS$PPH1]] <- results$summary["PP.H1.abf"]
	merged[[COLS$PPH2]] <- results$summary["PP.H2.abf"]
	merged[[COLS$PPH3]] <- results$summary["PP.H3.abf"]
	merged[[COLS$PPH4]] <- results$summary["PP.H4.abf"]
	merged[[COLS$PPH4.PPH3]] <- results$summary["PP.H4.abf"] / results$summary["PP.H3.abf"]
	merged[[COLS$METHOD]] <- "abf"

	############# Define Columns ##################
	core.columns <- c(COLS$CHR, COLS$F1, COLS$F2, COLS$PPH0, COLS$PPH1, COLS$PPH2, COLS$PPH3, COLS$PPH4, COLS$PPH4.PPH3, COLS$SNPID, COLS$SNP.PP.H4, COLS$Position)
	stats.columns <- c(COLS$P, COLS$BETA, COLS$SE, COLS$AF)

	############# Reorder Columns ##################
	merged <- merged[,  c(core.columns, paste0("F1_", stats.columns), paste0("F2_", stats.columns), COLS$METHOD)]
	asso <- merged[which.max(merged$SNP.PP.H4), c(core.columns, COLS$METHOD)]
	
	return( list(top=asso, full=merged) )
}