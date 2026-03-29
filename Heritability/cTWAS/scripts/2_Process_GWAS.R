################################################################################
################################################################################

#########################
######## UTILITY ########
#########################

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }

  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile)) {
    return(dirname(normalizePath(ofile)))
  }

  if (interactive()) {
    return(getwd())
  }

  stop("Cannot determine script directory")
}


################################################################################
################################################################################
check_col <- function(gwas) {
  nm <- names(gwas)

  if ("ID" %in% nm)   return("ID")
  if ("rsid" %in% nm) return("rsid")
  if ("SNPID" %in% nm) return("SNPID")
  # if ("id" %in% nm)   return("id")

  stop("No variant ID column found. Columns: ", paste(nm, collapse = ", "))
}

read_VCF <- function(path) {
  gwas <- VariantAnnotation::readVcf(path)
  gwas <- as.data.frame(gwasvcf::vcf_to_tibble(gwas))

  id_col <- check_col(gwas)
  

  z_snp <- ctwas::read_gwas(gwas,
                            id = id_col,
                            A1 = "ALT",
                            A2 = "REF",
                            beta = "ES",
                            se = "SE")

  return(z_snp)
}

################################################################################
################################################################################

# Z scores were reconstructed from GWAS summary statistics using reported effect
# sizes and two-sided p-values under a standard normal approximation.
# Specifically, Z was computed as sign(BETA) × Φ⁻¹(1 − P/2).
# This approach is equivalent to the Wald statistic (BETA/SE) in the limit
# of infinite precision.
# Minor numerical differences relative to SE-based Z scores are expected due
# to rounding of p-values and finite-precision arithmetic and do not affect
# downstream analyses.

check_GWAS_cols_for_z_score <- function(df) {
  required_cols <- c("BETA", "P")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Cannot compute Z score due to missing required columns: ",
      paste(missing_cols, collapse = ", "))
  }
}

calculate_z_score <- function(df) {
  # ------------------------------------------------------------------
  # Z-score reconstruction from GWAS summary statistics
  #
  # Some processed GWAS summary statistics do not include per-variant
  # standard errors (SE) or test statistics, but do report effect sizes
  # (BETA) and p-values (P). CTWAS requires per-variant Z scores rather
  # than raw effect sizes.
  #
  # Under standard GWAS practice (Wald test with two-sided p-values),
  # Z scores can be reconstructed as:
  #
  #   Z = sign(BETA) * qnorm(1 - P / 2)
  #
  # Assumptions:
  #   - P is a two-sided p-value
  #   - BETA corresponds to the effect allele (EA)
  #   - Test statistic is approximately normal under the null
  #
  # Minor numerical differences relative to SE-based Z scores are
  # expected due to p-value rounding and finite precision.
  # ------------------------------------------------------------------
  if (!"Z" %in% colnames(df)) {
    message("z score col not present in GWAS - attempting to recover z score")
    # stop("Z-score not present - file not GWAS summary statistic")
    # MODIFY DF IN PLACE
    check_GWAS_cols_for_z_score(df)
    # Replace zero p-values to avoid infinite z scores
    df[!is.na(P) & P == 0, P := 1e-300]
    # Compute z scores
    df[, Z := sign(BETA) * qnorm(1-P/2)]
    # Sanity checks (optional but recommended)
    #   cor(Z, BETA) should be positive
    #   Z should be approximately N(0,1) excluding tails
  } else {message("z score present in GWAS - proceeding")}
}

################################################################################
################################################################################
check_GWAS_cols_for_beta_score <- function(df) {

  cols <- names(df)
  has_Z_SE <- all(c("Z", "SE") %in% cols)
  has_Z_P  <- all(c("Z", "P") %in% cols)

  if (has_Z_SE) { return("Z_SE") }

  message("Z + SE not present - checking for Z + P")

  if (has_Z_P) {
    stop(
      "Z and P value are present. ",
      "SE could be estimated by providing GWAS_N. ",
      "This method is not yet implemented because GWAS_N is often unreliable."
    )
  }

  missing <- setdiff(c("Z", "SE", "P"), cols)

  stop(
    "Cannot compute BETA because required columns are missing: ",
    paste(missing, collapse = ", ")
  )
}

estimate_SE_from_N <- function(df, gwas_n) {
  df[, SE := 1 / sqrt(gwas_n)]
}

estimate_SE_from_N_and_MAF <- function(df, gwas_n) {
  df[, SE := 1 / sqrt(2 * gwas_n * MAF * (1 - MAF))]
}

beta_from_z_and_se <- function(df, gwas_n = NULL) {
  # gwas_n - dummy parameter
  df[, BETA := Z * SE]
  message("BETA added by Z * SE!")
  }

beta_from_Z_P <- function(df, gwas_n) {
  estimate_SE_from_N(df, gwas_n)
  message("SE calculated from P and gwas_n!")
  beta_from_z_and_se(df)
}

calculate_BETA_score <- function(df) {
  
  beta_calculator <- list(`Z_SE` = beta_from_z_and_se , `Z_P` = beta_from_Z_P)
  key <- check_GWAS_cols_for_beta_score(df)
  beta_calculator[key](df)
}
################################################################################
################################################################################

RSID_selector <- function(df) {
  # Trim spaces and convert empty strings to NA
  cols_to_clean <- c("RSID", "Original_RSID", "SNPID")
  # .SD = Subset of Data — all columns listed in .SDcols.
  # In this case, .SD will contain RSID, Original_RSID, and SNPID.
  # := is data.table syntax for by-reference assignment.
  # The parentheses around cols_to_clean indicate dynamic column names,
  # i.e., assign to the columns listed in that vector.
  # lapply(.SD, ...) applies the function to each column in
  df[, (cols_to_clean) := lapply(.SD, function(x) {
  # trims(x) -> removes leading and trailing whitespace (space, tabs)
  # from each string in the column.
  # fifelse(x == "", NA_character_, x) -> repleaces empty strings ("")
  # with NA (missing values).
  # fifelse is a fast, vectorised version of ifelse from data.table
  # NA_character_ ensures the column remains character type.
    x <- trimws(x)
    fifelse(x == "", NA_character_, x)
  # dplyr::coalesce(...) is a function that returns the first
  # non-missing value across its arguments.
  }), .SDcols = cols_to_clean]
  # Pick first non-missing ID
  df[, RSID_Sel := dplyr::coalesce(RSID, Original_RSID, SNPID)]
  return(df)
}

read_TSV <- function(path) {

  gwas <- data.table::fread(path, sep = "\t")

  # modify in place
  calculate_z_score(gwas)
  # calculate_BETA_score(gwas)
  RSID_selector(gwas)

  cols <- c("RSID_Sel", "NEA", "EA", "BETA", "Z")

  missing <- setdiff(cols, names(gwas))

  if (length(missing) > 0) {
    stop(
      "Columns missing from GWAS data - aborting. Missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  gwas <- gwas[, ..cols]

  setnames(
    gwas,
    old = cols,
    new = c("id", "A1", "A2", "ES", "z")
  )

  message("GWAS data read successfully")

  return(gwas)
}
################################################################################
################################################################################

# convert_to_ukb_varIDs <-function(varIDs, ref_format = "%s:%s_%s_%s") {
#     varID_list <- strsplit(varIDs, split = "_|:")
#     new_varIDs <- unlist(lapply(varID_list, function(x) {
#         if (length(x) >= 4) {
#             print(x)
#             sprintf(ref_format, x[1], x[2], x[3], x[4])
#         }
#         else {
#             x
#         }
#     }))
#     return(new_varIDs)
# }

check_z_score <- function(z_snp) {

  BAD_Z_SCORE <- !is.finite(z_snp$z)
  N_BAD_Z_SCORE <- sum(BAD_Z_SCORE)

  if (N_BAD_Z_SCORE > 0) {
    message("WARNING - dropping non-finite values from zsnp table")
    z_snp <- z_snp[!BAD_Z_SCORE]
    message("dropping ", N_BAD_Z_SCORE, " rows from gwas z scores")
  }

  return(z_snp)
}

convert_to_ukb_varIDs_fixed <- function(varIDs, ref_format = "%s:%s_%s_%s") {
  varID_list <- strsplit(varIDs, split = "_|:")
  vapply(seq_along(varID_list), function(i) {
    x <- varID_list[[i]]
    if (length(x) >= 4) sprintf(ref_format, x[1], x[2], x[3], x[4]) else varIDs[i]
  }, character(1))
}

prepare_gwas_znps <- function(GWAS_Path, SNP_map, key) {

  reader <- list(`VCF` = read_VCF, `TSV` = read_TSV)

  z_snp <- reader[[key]](GWAS_Path)
  message("Total SNPs in original file")
  message(capture.output(dim(z_snp)))

  z_snp <- check_z_score(z_snp)

  # Filter out missing, multiallelic, and strand ambiguous variants!!!
  # Convert variant IDs to match the refernce
  # Removes GWAS variants not present in LD reference!!!
  # Harmonises effect alleles with reference
  # Returns a read-to-use GWAS z-score table for cTWAS
  z_snp <- ctwas::preprocess_z_snp(z_snp, SNP_map,
                                   drop_multiallelic = TRUE,
                                   drop_strand_ambig = TRUE,
                                   varID_converter_fun = NULL)
                                  #  varID_converter_fun = convert_to_ukb_varIDs_fixed)

  message("Total SNPs in harmonised file")
  message(capture.output(dim(z_snp)))
  invisible(z_snp)
}

################################################################################
################################################################################

extract_name <- function(GWAS_Path) {

  parent_dir <- basename(dirname(GWAS_Path))
  print("parent_dir")
  print(parent_dir)
  study_name <- gsub("\\.(db|csv|tsv|vcf|gz)", "", basename(GWAS_Path))
  print("study_name")
  print(study_name)
  label <- paste(parent_dir, study_name, sep = "_")
  print("label")
  print(label)

  return(list(`parent_dir` = parent_dir,
              `study_name` = study_name,
              `label` = label))

}

################################################################################
################################################################################



################################################################################
################################################################################

extract_file_type_from_inpath <- function(inputpath) {
  compression_ext <- c("gz", "bgz")

  fname <- basename(inputpath)

  # Extract last extension
  ext <- tolower(tools::file_ext(fname))

  # if compressed, strip compression and re-check extension
  if (ext %in% compression_ext) {
    fname2 <- sub(paste0("\\.", ext, "$"), "", fname, ignore.case=TRUE)
    ext <- tolower(tools::file_ext(fname2))
  }

  ftype_key_dispensor <- list(
    "vcf" = "VCF",
    "tsv" = "TSV",
    "txt" = "TSV",
    "csv" = "CSV"
  )

  key <- ftype_key_dispensor[[ext]]

  stop_txt <- paste0("Unsupported file type: ", fname,
                     " (detected ext='", ext, "'). ",
                     "Allowed: .vcf, .tsv, .txt, .csv ",
                     "optional .gz/.bgz. compressions")

  if (is.null(key)) { stop(stop_txt) }
  print(paste("Reading file of type", ext, "file name - ", fname))
  return(key)
}

################################################################################
################################################################################
################################################################################
################################################################################

main <- function(GWAS_Path, SNP_map_path, genome_version, outpath, disease) {

  ftype <- extract_file_type_from_inpath(GWAS_Path)

  print("Executing Script")
  z_snp <- prepare_gwas_znps(GWAS_Path, readRDS(SNP_map_path), ftype)

  print("Extracting name")
  disease <- extract_name(GWAS_Path)[['parent_dir']]
  study_name <- extract_name(GWAS_Path)[['study_name']]

  print("Saving file")
  save_RDS(z_snp, outpath)

}

if (sys.nframe() == 0) {
  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("ctwas")
    library("argparse")
    library("data.table")
  })
  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--GWAS_Path", help = "region path", required = TRUE)

  parser$add_argument("--SNP_map_path", help = "SNP map for reference", required = TRUE)

  parser$add_argument("--genome_version", help = "genome_version",
                      required = FALSE,
                      default = "b38")

  parser$add_argument("--chromosome", type = "integer",
                      help = "chromosme, for subsetting", required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath", help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--disease", help = "file type to read",
                      required = FALSE,
                      default = NULL)

  ######################################################

  args <- parser$parse_args()

  ######################################################
  main(GWAS_Path = args$GWAS_Path,
       SNP_map_path = args$SNP_map,
       genome_version = args$genome_version,
       outpath = args$outpath,
       disease = args$disease)
}