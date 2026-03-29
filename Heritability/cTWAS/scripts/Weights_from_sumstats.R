################################################################################
################################################################################

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

configure <- function(SNPID_FORMAT, CONVERT_ID_FORMAT, DB_PATH) {

  if (CONVERT_ID_FORMAT && !file.exists(DB_PATH)) (
    stop("WARNING - Selected to convert SNPID to RSID but DB_PATH does not exist.")
  )

  if (SNPID_FORMAT == "SNPID" && CONVERT_ID_FORMAT) { CONVERT_ID_FORMAT_COL = "RSID" }

  if (SNPID_FORMAT == "RSID" && CONVERT_ID_FORMAT) { CONVERT_ID_FORMAT_COL = "RSID" }

  return( list(
               `CONVERT_ID_FORMAT` = CONVERT_ID_FORMAT,
               `SNPID_FORMAT` = SNPID_FORMAT,
               `CONVERT_ID_FORMAT_COL` = CONVERT_ID_FORMAT_COL
               )
               )
}

################################################################################
################################################################################


read_gtex_eqtls <- function(path) {
  DT <- fread(path)

  extract_cols <- c("GENEID", "rsid", "SNPID", "ref", "alt", "BETA", "P")

  DT <- DT[, ..extract_cols]

  colnames(DT) <-  c("gene",
                     "rsid",
                     "varID",
                     "ref_allele",
                     "eff_allele",
                     "weight",
                     "pval")

  # Add chromosome column
  DT[, Chromosome := paste0("chr", sub(":.*", "", DT$varID))]
  return(DT)
}

################################################################################
################################################################################

read_dice_eqtls <- function(path) {
  # # # # # # # # # # # # # # # # # # 
  DT <- fread(path)
  # # # # # # # # # # # # # # # # # # 
  # TO REPLACE for DICE CELLTYPES!!
  if ("rsid" %in% colnames(DT) || "RSID" %in% colnames(DT)) {
  expected_cols <- c("GENEID", "rsid", "SNPID", "ref", "alt", "BETA", "P")
  new_names <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight", "pval")
  } else { 
  expected_cols <- c("GENEID", "SNPID", "ref", "alt", "BETA", "P")
  new_names <- c("gene", "varID", "ref_allele", "eff_allele", "weight", "pval")
  }
  # # # # # # # # # # # # # # # # # # 
  DT <- DT[, ..expected_cols]
  # TO REPLACE for DICE CELLTYPES!!
  colnames(DT) <-  new_names
  # # # # # # # # # # # # # # # # # # 
  DT[, Chromosome := paste0("chr", sub(":.*", "", DT$varID))]
  return(DT)
  # # # # # # # # # # # # # # # # # # 
}

################################################################################
################################################################################

prepare_gene_table <- function(gtf, gene_list) {

  df <- S4Vectors::mcols(gtf)[, c("gene_id", "gene_name", "gene_type")]
  df <- df[df$gene_id %in% gene_list, ]
  df <- unique(df)
  df <- as.data.table(as.data.frame(df))

  setnames(
           df,
           old = c("gene_id", "gene_name", "gene_type"),
           new = c("gene", "genename", "gene_type")
           )

  return(df)
}

################################################################################
################################################################################

assemble_chromosome_ld_path_map <- function(low_level_chrom_path_hblock) {
  file_paths <- list.files(low_level_chrom_path_hblock, full.names=TRUE)

  ld_idx <- grepl("\\.RDS$", file_paths)
  snp_idx <- grepl("\\.Rvar$", file_paths)

  ld_files <- file_paths[ld_idx]
  snp_files <- file_paths[snp_idx]

  ld_names <- basename(ld_files)
  hblocks <- sub(".*\\.(.*)\\..*$", "\\1", ld_names)
  chroms <- sapply(strsplit(ld_names, "\\.|_"), `[`, 5)

  # sanity check
  if(length(ld_files) != length(snp_files)) {
    warning("Number of LD files does not match number of SNP files")
  }
  # Build data.table
  dt <- data.table::data.table(
    Chromosome = chroms,
    region_id = paste0(chroms, "_", hblocks),
    LD_file = ld_files,
    SNP_file = snp_files
  )
  return(dt)
}

build_path_frame_chrom_all <- function(LD_directory) {
  print(LD_directory)
  chrom_dirs <- list.dirs(LD_directory, recursive = FALSE, full.names = TRUE)
  dt_list <- lapply(chrom_dirs, assemble_chromosome_ld_path_map)
  dt_list <- rbindlist(dt_list)
  dt_by_chromosome <- split(dt_list, by = "Chromosome")
  return(dt_by_chromosome)
}

################################################################################
################################################################################

extract_id_from_db <- function(DB_PATH, SNPIDs, col = "SNPID", CONVERT_ID_FORMAT_COL) {

  con <- DBI::dbConnect(RSQLite::SQLite(), DB_PATH)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  if (!col %in% DBI::dbListFields(con, "Variants")) {
    stop("Error: column '", col, "' not found in Variants table")
  }

  if (length(SNPIDs) == 0) {
  stop("Error: no SNPIDs have been passed into 'extract_id_from_db' function.")
  }

  snps.use <- paste(DBI::dbQuoteString(con, SNPIDs), collapse = ",")

  query <- glue::glue(
    "SELECT * FROM Variants WHERE {col} IN ({snps.use})"
  )

  dbsnp <- data.table(data.frame(DBI::dbGetQuery(con, query)))
  if (CONVERT_ID_FORMAT_COL == "RSID" && CONVERT_ID_FORMAT_COL %in% colnames(dbsnp)) {
    setnames(dbsnp, "RSID", "rsid")
    CONVERT_ID_FORMAT_COL <- "rsid"
  }
  # if ("RSID" %in% colnames(dbsnp)) {
  #   data.table::setnames(dbsnp, "RSID", "rsnumber_156")
  # }
  print(dbsnp)
  cols <- c(col, CONVERT_ID_FORMAT_COL)
  return(
         list(dbsnp[, ..cols], 
              CONVERT_ID_FORMAT_COL)
              )
}

merge_qtl_table_and_db_ <- function(DB_PATH, qtl_table, SNPID_FORMAT, CONVERT_ID_FORMAT_COL) {

  dbsnp_l <- extract_id_from_db(
    DB_PATH,
    qtl_table[["varID"]],
    col = SNPID_FORMAT,
    CONVERT_ID_FORMAT_COL
  )

  dbsnp <- dbsnp_l[[1]]
  CONVERT_ID_FORMAT_COL <- dbsnp_l[[2]]

  dbsnp <- merge(
    qtl_table,
    dbsnp,
    by.x = "varID",
    by.y = SNPID_FORMAT
  )

  setnames(
    dbsnp,
    old = c("varID", CONVERT_ID_FORMAT_COL),
    new = c(SNPID_FORMAT, "varID")
  )

  return(dbsnp)
}

extract_LD_info_from <- function(matrix_path, var_path, snps.keep) {
  M <- readRDS(matrix_path)
  var <- read.table(var_path, header = TRUE, stringsAsFactors = TRUE)

  stopifnot(nrow(M) == ncol(M))

  keep <- match(snps.keep, var$id)
  keep <- keep[!is.na(keep)]

  if (!length(keep)) {
    return(data.table(RSID1 = character(),
                      RSID2 = character(),
                      LD = numeric()))
  }

  LDsub <- M[keep, keep, drop = FALSE]
  snp_subset <- var$id[keep]

  ind <- which(upper.tri(LDsub, diag = TRUE), arr.ind = TRUE)
  data.table(
    RSID1 = snp_subset[ind[, 1]],
    RSID2 = snp_subset[ind[, 2]],
    LD    = LDsub[ind]
  )
}

################################################################################
################################################################################

load_table <- function(qtl_path, qtl_data_base) {
  # Load QTL table
  read_selector <- list("GTEX" = read_gtex_eqtls, "DICE" = read_dice_eqtls)

  qtl_table <- read_selector[[qtl_data_base]](qtl_path)
  qtl_table[, Chromosome := paste0("chr", sub(":.*", "", qtl_table$varID))]
  return(qtl_table)
}

extract_LD_info_for_singular_chrom <- function(chr, LD_path_map_by_chr, snps_by_chr) {

  message("Processing ", chr)
  ld_dt <- LD_path_map_by_chr[[chr]]
  snps_chr <- snps_by_chr[[chr]]

  tables <- lapply(
    seq_len(nrow(ld_dt)),
        function(i) {
          extract_LD_info_from(
            ld_dt$LD_file[i],
            ld_dt$SNP_file[i],
            snps_chr
          )
        }
      )

  return(data.table::rbindlist(tables))
}

qtls_by_chromosome <- function(qtl_table, LD_path_map_by_chr) {
  # Extract LD matrix paths

  qtl_table <- qtl_table[Chromosome %chin% names(LD_path_map_by_chr)]

  qtl_table_by_chr <- qtl_table[
    , .(varID = list(unique(varID))), by = Chromosome ]

  snps_by_chr <- setNames(
    qtl_table_by_chr$varID,
    qtl_table_by_chr$Chromosome
  )

  return(snps_by_chr)

}

create_predict_df_from_QTLs_wrapper <- function(qtl_table,
                                                gene_table,
                                                cov_table,
                                                outdir,
                                                qtl_path) {

  outname <- tools::file_path_sans_ext(basename(qtl_path))
  # Build PredictDB
  ctwas::create_predictdb_from_QTLs(
    weight_table = qtl_table,
    gene_table   = gene_table,
    cov_table    = cov_table,
    use_top_QTL  = FALSE,
    select_by    = "weight",
    outputdir    = outdir,
    outname      = outname
  )
  return(outname)
}

extract_LD_info_mc_main <- function(LD_directory, qtl_table) {
  # Extract LD matrix paths
  # Split LD and snps by chromosome
  LD_path_map_by_chr <- build_path_frame_chrom_all(LD_directory)
  snps_by_chr <- qtls_by_chromosome(qtl_table, LD_path_map_by_chr)

  # Prepare LD table
  cov_table <- mclapply(
    names(LD_path_map_by_chr),

  extract_LD_info_for_singular_chrom,
  LD_path_map_by_chr = LD_path_map_by_chr,
  snps_by_chr = snps_by_chr,
  mc.cores = 10
  )

  return(data.table::rbindlist(cov_table))
}

main <- function(qtl_path,
                 qtl_data_base,
                 gtfpath,
                 LD_directory,
                 outdir,
                 SNPID_FORMAT,
                 CONVERT_ID_FORMAT_COL,
                 CONVERT_ID_FORMAT,
                 DB_PATH) {

  # Load QTL table
  qtl_table <- load_table(qtl_path, qtl_data_base)

  if (CONVERT_ID_FORMAT) {
  qtl_table <- merge_qtl_table_and_db_(DB_PATH,
                                       qtl_table,
                                       SNPID_FORMAT,
                                       CONVERT_ID_FORMAT_COL)
  }

  # Prepare gene table using imported GTF
  gtf <- rtracklayer::import(gtfpath, feature.type = "gene")
  gene_table <- prepare_gene_table(gtf, unique(qtl_table$gene))

  # Prepare LD table
  cov_table <- extract_LD_info_mc_main(LD_directory = LD_directory,
                                       qtl_table = qtl_table)
  # Add gene information to LD table
  cov_table <- merge(cov_table,
                     qtl_table[, c("gene", "rsid")],
                     by.x = "RSID1",
                     by.y = "varID",
                     allow.cartesian = TRUE)
  # Rename relevant columns post merging tables
  data.table::setnames(cov_table,
                       old = c("gene", "LD"),
                       new = c("GENE", "VALUE"))

  # Generate database information
  outname <- create_predict_df_from_QTLs_wrapper(qtl_table,
                                                 gene_table,
                                                 cov_table,
                                                 outdir, qtl_path)

  # Return path of created DB
  file.path(outdir, paste0(outname, ".db"))

}

################################################################################
################################################################################

# main <- function(qtl_data_base, qtl_paths, gtfpath, outdir, ncores = 1) {

#   # Import GTF once
#   gtf <- rtracklayer::import(gtfpath, feature.type = "gene")

#   qtl_paths <- path_parser(qtl_paths, qtl_data_base)
#   print(qtl_paths)
  
#   out <- parallele_processes(qtl_paths, qtl_data_base, gtf, outdir, ncores)
#   return(out)
# }

################################################################################
################################################################################

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("DBI")
    library("glue")
    library("ctwas")
    library("RSQLite")
    library("parallel")
    library("argparse")
    library("data.table")
    library("rtracklayer")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--qtl_data_base",
                      help = "qtl database origin - e.g. DICE or GTEX",
                      required = TRUE,
                      choices = c("DICE", "GTEX"))

  parser$add_argument("--qtl_path",
                      nargs = "*",
                      help = "path to qtl yaml file or n number of paths",
                      required = TRUE)

  parser$add_argument("--gtfpath",
                      help = "path to GTF file.",
                      required = TRUE)

  parser$add_argument("--LD_directory",
                      help = "path to LD files for all chroms",
                      required = TRUE)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--SNPID_FORMAT",
                      required = FALSE,
                      choices = c("SNPID", "RSID"),
                      default = "SNPID")

  parser$add_argument("--CONVERT_ID_FORMAT_COL",
                      required = FALSE,
                      choices = c("SNPID", "RSID"),
                      default = "SNPID")

  parser$add_argument("--CONVERT_ID_FORMAT",
                      required = FALSE,
                      default = FALSE,
                      action = 'store_true')

  parser$add_argument("--DB_PATH",
                      required = FALSE)

  # parser$add_argument("--genome_version",
  #                     help = "genome_version",
  #                     required = TRUE)

  ######################################################
  args <- parser$parse_args()
  ######################################################

  config <- configure(
    args$SNPID_FORMAT,
    args$CONVERT_ID_FORMAT,
    args$DB_PATH
    )

  main(qtl_path = args$qtl_path,
       qtl_data_base = args$qtl_data_base,
       gtfpath = args$gtfpath,
       LD_directory = args$LD_directory,
       outdir = args$outpath,
       SNPID_FORMAT = config$SNPID_FORMAT,
       CONVERT_ID_FORMAT_COL = config$CONVERT_ID_FORMAT_COL,
       CONVERT_ID_FORMAT = args$CONVERT_ID_FORMAT,
       DB_PATH = args$DB_PATH)
}

# DB_PATH <- "/home/jrocha/BioAdHoc/reference_panels/dbSNP/hg38/dbSNP156.db"
################################################################################
################################################################################
# path <- "/home/jottensmeier/BioAdHoc/Projects/CTCF/Analysis_GRCh38_New/Intermediary/Input_Clean/filt_eQTL_catalogue/filtered/eQTL_Catalogue_1e_neg4_GTEx_ge_adipose_subcutaneous_cleaned.tsv"
# outdir <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/c_TWAS/results/eqtl_weights"
# tazle <- read_gtex_eqtls(path)
# gene_table <- prepare_gene_table(gtfp, unique(table$gene))
# gtfp <- paste0("/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/",
#                "reference/GRCh38-2020-A_build/gtf/",
#                "gencode.v32.primary_assembly.annotation.filtered.gtf")
# yaml::read_yaml()
################################################################################
################################################################################
# qtlPoc <- "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/c_TWAS/results/eqtl_weights/eQTL_Catalogue_1e_neg4_GTEx_ge_adipose_subcutaneous_cleaned.db"
# con <- dbConnect(RSQLite::SQLite(), dbname = qtlPoc)
# table <- dbReadTable(con, "weights")
# tab2 <-  dbReadTable(con, "extra")
# dbDisconnect(con)
################################################################################
################################################################################
# rename_cols <- function(table, mapper) {
#   cols <- colnames(table)
#   to_replace <- cols %in% names(mapper)
#   cols[to_replace] <- unname(mapper[cols[to_replace]])
#   colnames(table) <- cols
#   invisible(table)
# }
# ################################################################################
# ################################################################################
# read_gtex_eqtls <- function(path) {
#   table <- read.csv(path, sep = "\t")
#   mapper <- c("GENEID" = "gene",
#               "rsid" = "rsid",
#               "SNPID" = "varID",
#               "ref" = "ref_allele",
#               "alt" = "eff_allele",
#               "BETA" = "weight",
#               "P" = "pval")
#   table <- rename_cols(table, mapper)
#   table <- table[, unname(mapper), drop = FALSE]
#   invisible(table)
# }





################################################################################
################################################################################

#' @title Computes LD for weight variants using reference LD
#'
#' @param weights a list of preprocessed weights.
#'
#' @param region_info a data frame of region definitions.
#'
#' @param LD_map a data frame with filenames of LD matrices and SNP information for the regions.
#' Required when \code{load_predictdb_LD = FALSE}.
#'
#' @param snp_map a list of SNP-to-region map for the reference.
#' If NUll, it will reads SNP info from the "SNP_file" column of LD_map.
#'
#' @param LD_format file format for LD matrix. If "custom", use a user defined
#' \code{LD_loader_fun()} function to load LD matrix.
#'
#' @param LD_loader_fun a user defined function to load LD matrix when \code{LD_format = "custom"}.
#'
#' @param snpinfo_loader_fun a user defined function to load SNP information file,
#' if SNP information files are not in standard cTWAS reference format.
#'
#' @param ncore The number of cores used to parallelize computation.
#'
#' @importFrom parallel mclapply
#' @importFrom Matrix bdiag
#' @importFrom logging loginfo
#'
#' @return a list of processed weights, with LD of weight variants included.
#'
#' @export
# compute_weight_LD_from_ref <- function(weights,
#                                        region_info,
#                                        LD_map,
#                                        snp_map = NULL,
#                                        LD_format = c("rds", "rdata", "mtx", "csv", "txt", "custom"),
#                                        LD_loader_fun = NULL,
#                                        snpinfo_loader_fun = NULL,
#                                        ncore = 1) {

#   LD_format <- match.arg(LD_format)

#   if (!inherits(weights,"list"))
#     stop("'weights' should be a list!")

#   if (!inherits(LD_map,"data.frame"))
#     stop("'LD_map' should be a data frame!")

#   weight_info <- lapply(names(weights), function(x){
#     as.data.frame(weights[[x]][c("chrom", "p0","p1", "molecular_id", "weight_name", "type","context")])})
#   weight_info <- do.call(rbind, weight_info)
#   weight_info$weight_id <- paste0(weight_info$molecular_id, "|", weight_info$weight_name)
#   # get the regions overlapping with each gene
#   for (k in 1:nrow(weight_info)) {
#     chrom <- weight_info[k, "chrom"]
#     p0 <- weight_info[k, "p0"]
#     p1 <- weight_info[k, "p1"]
#     idx <- which(region_info$chrom == chrom & region_info$start <= p1 & region_info$stop > p0)
#     weight_info[k, "region_id"] <- paste(sort(region_info[idx, "region_id"]), collapse = ",")
#   }

#   # compute LD for weight variants on each chromosome
#   chrs <- sort(unique(weight_info$chrom))
#   for (b in chrs) {
#     loginfo("Computing LD for variants in weights on chr%s", b)
#     weightinfo <- weight_info[weight_info$chrom == b, ]
#     if (nrow(weightinfo) > 0) {
#       weight_region_ids <- names(sort(-table(weightinfo$region_id)))
#       weight_LD_list <- mclapply_check(weight_region_ids, function(x){
#         # load the R_snp and SNP info for the region
#         # and extract LD for the weight variants
#         curr_region_LD_list <- list()
#         curr_region_ids <- unlist(strsplit(x, ","))
#         curr_region_idx <- match(curr_region_ids, LD_map$region_id)

#         LD_matrix_files <- unlist(strsplit(LD_map$LD_file[curr_region_idx], split = ","))
#         stopifnot(all(file.exists(LD_matrix_files)))

#         if (length(LD_matrix_files) > 1) {
#           R_snp <- lapply(LD_matrix_files, load_LD, format = LD_format, LD_loader_fun = LD_loader_fun)
#           R_snp <- suppressWarnings(as.matrix(bdiag(R_snp)))
#         } else {
#           R_snp <- load_LD(LD_matrix_files, format = LD_format, LD_loader_fun = LD_loader_fun)
#         }

#         # load SNP info of the region
#         # if snp_map is available, reads SNP info from snp_map;
#         # otherwise, reads SNP info from the "SNP_file" column of LD_map.
#         if (!is.null(snp_map)){
#           snpinfo <- as.data.frame(rbindlist(snp_map[curr_region_ids], idcol = "region_id"))
#         } else {
#           SNP_info_files <- LD_map$SNP_file[curr_region_idx]
#           stopifnot(all(file.exists(SNP_info_files)))
#           snpinfo <- read_snp_info_files(SNP_info_files, snpinfo_loader_fun = snpinfo_loader_fun)
#         }

#         rownames(R_snp) <- snpinfo$id
#         colnames(R_snp) <- snpinfo$id
#         weight_ids <- weightinfo[weightinfo$region_id == x, "weight_id"]

#         for (weight_id in weight_ids) {
#           wgt_snp_ids <- rownames(weights[[weight_id]]$wgt)
#           R_wgt <- R_snp[wgt_snp_ids, wgt_snp_ids, drop=FALSE]
#           curr_region_LD_list[[weight_id]] <- R_wgt
#         }
#         curr_region_LD_list
#       }, mc.cores = ncore, stop_if_missing = TRUE)

#       weight_LD_list <- unlist(weight_LD_list, recursive = FALSE)
#       for(weight_id in names(weight_LD_list)){
#         weights[[weight_id]][["R_wgt"]] <- weight_LD_list[[weight_id]]
#       }
#     }
#   }
#   return(weights)
# }
