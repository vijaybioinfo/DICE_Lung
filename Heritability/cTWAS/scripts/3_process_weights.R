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

#########################
######## UTILITY ########
#########################

#########################
##### CONFIG PARSER #####
#########################

read_config <- function(path, key = NULL) {
  config <- yaml::read_yaml(path)
  if (!is.null(key) && key %in% names(config)) {
    config <- config[[key]]
  }
  return(config)
}

load_weight_spec <- function(weight_file,
                             model = NULL,
                             type = NULL,
                             context = NULL,
                             key = "prediction_models") {

  is_yaml <- is.character(weight_file) &&
    tools::file_ext(weight_file) %in% c("yaml", "yml")

  is_config_mode <- is_yaml && is.null(type) && is.null(context)
  is_inline_mode <- !is.null(type) && !is.null(context)

  if (is_config_mode == is_inline_mode) {
    stop("Specify either YAML config OR tyep/context (exclusively)")
  }

  if (is_config_mode) {
    print('Running Config Mode for weight harmonisation')
    config <-  read_config(weight_file, key)

  } else {
    print('Running Inline Mode for weight harmonisation')
    config <- list(list(
    model = model,
    type = type,
    context = context,
    path = weight_file
  ))

  return(config)
  }
}

#########################
######## WRITER  ########
#########################

# refer to UTILS

################################################################################
################################################################################

#########################
######## WRAPPER ########
#########################

process_weights_wrapper <- function(weight_file,
                                    region_info,
                                    snp_map,
                                    gwas_snp_ids,
                                    type,
                                    context,
                                    LD_map = NULL,
                                    load_predictdb_LD = TRUE) {

  weights <- ctwas::preprocess_weights(
    weight_file = weight_file,
    region_info = region_info,
    gwas_snp_ids = gwas_snp_ids,
    snp_map = snp_map,
    type = type,
    context = context,
    weight_name = paste0(type, "_", context),
    weight_format = "PredictDB",
    drop_strand_ambig = TRUE,
    filter_protein_coding_genes = FALSE,
    scale_predictdb_weights = FALSE,
    ##################################
    ## !!!Check these parameters!!! ##
    ##################################
    LD_map = LD_map,
    load_predictdb_LD = load_predictdb_LD,
    include_weight_LD = TRUE,
    ncore = 2,
    LD_format = "rds"
    # top_n_snps = 500
    # weight_format = "PredictDB"
  )
  invisible(weights)
}

################################################################################
################################################################################

#########################
####### Processing ######
#########################

processing <- function(weight_file,
                       region_info,
                       snp_map,
                       gwas_snp_ids,
                       prediction_model,
                       type,
                       context,
                       outpath,
                       LD_map_path = NULL) {

  # Determine LD source
  if (file.exists(outpath) && file.info(outpath)$size > 0) {
    message("Harmonised weights file exists - skipping time intensive harmonisation step.")
    message(outpath)
    return(readRDS(outpath))
    }

  LD_map <- if (!is.null(LD_map_path)) readRDS(LD_map_path) else NULL

  print(weight_file)
  print(glue('Null status of LD map path {is.null(LD_map_path)}.'))
  print(glue('Null status of LD map {is.null(LD_map)}.'))

  ## Calculate weights
  weight <- process_weights_wrapper(
    weight_file = weight_file,
    region_info = region_info,
    snp_map = snp_map,
    gwas_snp_ids = gwas_snp_ids,
    type = paste(prediction_model, type, sep = "_"),
    context = context,
    LD_map = LD_map,
    load_predictdb_LD = is.null(LD_map)
  )

  print("Processing finished")
  ## Save weights
  # label <- paste(prediction_model, type, context, sep = "_")
  save_RDS(weight, outpath)
  invisible(weight)
}

################################################################################
################################################################################

#########################
######## UTILITY ########
#########################

make_worker <- function(region_info,
                        snp_map,
                        gwas_snp_ids,
                        outpath,
                        LD_map_path = NULL) {

  LD_map <- if (!is.null(LD_map_path)) readRDS(LD_map_path) else NULL

  function(weight_file, context, type, prediction_model = NULL) {
    processing(
      weight_file = weight_file,
      region_info = region_info,
      snp_map = snp_map,
      gwas_snp_ids = gwas_snp_ids,
      prediction_model = prediction_model,
      type = type,
      context = context,
      outpath = outpath,
      LD_map_path = LD_map_path
    )
  }
}

#########################
#######  PARALLEL #######
#########################

single <- function(worker, weights, ncores) {
  print(weights)
  worker(
    weight_file = weights[[1]]$path,
    context = weights[[1]]$context,
    type = weights[[1]]$type,
    prediction_model = weights[[1]]$model
  )
}

parallel <- function(worker, weights, ncores) {
  # Multiple weight files.
  print("##################################")
  print(weights)
  print("##################################")
  multigroup_weights <- parallel::mclapply(
    weights, function(args) {
      print("##################################")
      print(args)
      print("##################################")
      worker(weight_file = args$path,
             context = args$context,
             type = args$type,
             prediction_model = args$model)
             },
    mc.cores = ncores
  )
  return(multigroup_weights)
}

dispatcher <- list(`single` = single,
                   `multiprocessing` = parallel)

################################################################################
################################################################################

main <- function(weight_files,
                 region_path,
                 chromosome,
                 snp_map_path,
                 gwas_snpid_path,
                 LD_map_path = NULL,
                 outpath,
                 ncores = 1)  {

  # Output directory.
  # study_spec_outpath <- generate_outpath(
  #   outpath,
  #   file.path('harmonised_weights',
  #             extract_base_dir_label(gwas_snpid_path))
  # )

  # Load GWAS SNP IDs.
  gwas_snp_ids <- readRDS(gwas_snpid_path)$id

  # Prepare regions.
  region_info <- prepare_region_info(region_path,
                                     chromosome)

  # Load SNP map.
  snp_map <- readRDS(snp_map_path)

  print(glue("Harmonising {length(weight_files)} weight files."))

  # Create worker (closure)
  worker <- make_worker(region_info,
                        snp_map,
                        gwas_snp_ids,
                        outpath,
                        LD_map_path)

  print("Worker Loaded - beginning weight harmonisation")
  # Run in parallel
  proc_type <- if (length(weight_files) == 1) "single" else "multiprocessing"
  weights <- dispatcher[[proc_type]](worker, weight_files, ncores)

}

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("ctwas")
    library("argparse")
    library("parallel")
    library("data.table")
  })

  ######################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--weight_files",
                      help = "weight file path per condition",
                      required = TRUE)

  parser$add_argument("--region_path",
                      help = "LD dir",
                      required = TRUE)

  parser$add_argument("--chromosome",
                      type = "integer",
                      help = "chromosme, for subsetting",
                      required = FALSE, default = NULL)

  parser$add_argument("--snp_map_path",
                      help = "snp_map",
                      required = TRUE)

  parser$add_argument("--gwas_snpid_path",
                      help = "GWAS SNPID references snp ids",
                      required = TRUE)

  parser$add_argument("--prediction_model",
                      help = "tissue type",
                      default = NULL,
                      required = FALSE)

  parser$add_argument("--type",
                      help = "qtl type",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--context",
                      help = "tissue type",
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--LD_map_path",
                      help = paste0("A data frame with filenames of",
                                    "LD matrices and SNP information",
                                    "for the regions. Required when",
                                    "load_predictdb_LD = FALSE"),
                      required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  parser$add_argument("--ncores", help = "number of processing cores",
                      default = 1,
                      type = "integer",
                      required = FALSE)

  ######################################################
  args <- parser$parse_args()
  ######################################################
  print(args$weight_files)

  # Parse running mode - config or inline.
  weight_files <- load_weight_spec(
    weight_file = args$weight_files,
    model = args$prediction_model,
    type = args$type,
    context = args$context,
    key = "prediction_models"
  )

  main(
    weight_file = weight_files,
    region_path = args$region_path,
    chromosome = args$chromosome,
    snp_map_path = args$snp_map_path,
    gwas_snpid_path = args$gwas_snpid_path,
    LD_map_path = args$LD_map_path,
    outpath = args$outpath,
    ncores = args$ncores
  )
  ######################################################
}

# weight_files <- load_weight_spec(
#     weight_files,
#     prediction_model,
#     type,
#     context,
#     key = "prediction_models"
#   )

# # weights <- weight_files

  #   weight_files <- load_weight_spec(
  #   weight_file = weight_files,
  #   model = prediction_model,
  #   type = type,
  #   context = context,
  #   key = "prediction_models"
  # )

#  main(
#     weight_file = weight_files,
#     region_path = region_path,
#     chromosome = NULL,
#     snp_map_path = snp_map_path,
#     gwas_snpid_path = gwas_snpid_path,
#     LD_map_path = LD_map_path,
#     outpath = outpath,
#     ncores = ncores
#   )

#   # for testing
#   zsnpid <- read.csv("/home/jottensmeier/BioAdHoc/Projects/DICE_LUNG/c_TWAS/results/harmonised_gwas_zsnps/LDL_Top_GTEX_eqtl_ukb-d-30780_irnt.vcf.gz.rds")$id
#   region_path <- "/mnt/biohome/jottensmeier/miniforge3/envs/cTWAS/lib/R/library/ctwas/extdata/ldetect/EUR.b38.ldetect.regions.RDS"
#   chromosome <- 16
#   region_info <- prepare_region_info(region_path, chromosome)
#   snp_map <- "/home/jottensmeier/BioAdHoc/Projects/DICE_LUNG/c_TWAS/results/snp_map/ukb_b38_0.1_chr16.R_snp.rds"
#   snp_map <- readRDS(snp_map)

#   weight_file <- system.file("extdata/sample_data", "expression_Liver.db", package = "ctwas")
# #   weight_adipose_file <- system.file("extdata/sample_data", "expression_Adipose_Subcutaneous.db", package = "ctwas")
#   context <- "Liver"
#   type <- "eqtl"

#   label = paste(context, type, sep = "_")

#   weights <- process_weights_wrapper(weight_file, region_info, zsnpid,
#                                      snp_map, context, type)




# weight_name = paste0(context, "_", type)
# weight_format = "PredictDB"
# drop_strand_ambig = TRUE
# filter_protein_coding_genes = TRUE
# scale_predictdb_weights = TRUE
# load_predictdb_LD = TRUE
# include_weight_LD = TRUE
# fusion_method = "enet"
# fusion_genome_version = NA
# top_n_snps = NULL
# LD_format = "rds"
# LD_loader_fun = NULL
# snpinfo_loader_fun = NULL
# ncore = 1
# logfile = NULL
# verbose = TRUE


#   weights <- ctwas:::mclapply_check(weight_molecular_ids, function(molecular_id){
#     ctwas:::process_weights(molecular_id,
#                     type = type,
#                     context = context,
#                     weight_name = weight_name,
#                     weight_table = weight_table,
#                     cov_table = cov_table,
#                     snp_info = snp_info,
#                     weight_format = weight_format,
#                     top_n_snps = top_n_snps,
#                     drop_strand_ambig = drop_strand_ambig,
#                     scale_predictdb_weights = scale_predictdb_weights,
#                     verbose = verbose)
#   }, mc.cores = ncore)


# process_weights <- function (molecular_id,
#                              type,
#                              context,
#                              weight_name,
#                              weight_table, 
#                              cov_table,
#                              snp_info,
#                              weight_format = c("PredictDB", "FUSION"),
#                              top_n_snps = NULL,
#                              drop_strand_ambig = TRUE,
#                              scale_predictdb_weights = TRUE, 
#                              verbose = FALSE)
# {
#     if (weight_format == "FUSION") {
#         scale_predictdb_weights <- FALSE
#     }
#     target_header <- c("chrom", "id", "pos", "alt", "ref")
#     if (!all(target_header %in% colnames(snp_info))) {
#         stop("snp_info needs to contain the following columns: ",
#             paste(target_header, collapse = " "))
#     }
#     if (verbose)
#         loginfo("Processing weight for %s...", molecular_id)
#     g.weight_table <- weight_table[weight_table$gene == molecular_id,
#         ]
#     wgt.matrix <- as.matrix(g.weight_table[, "weight", drop = FALSE])
#     rownames(wgt.matrix) <- g.weight_table$rsid
#     chrom <- sapply(strsplit(g.weight_table$varID, "_"), "[[", 
#         1)
#     chrom <- unique(parse_number(chrom))
#     if (length(chrom) > 1) {
#         stop(sprintf("More than one 'chrom' in weight of %s!", 
#             molecular_id))
#     }
#     snp_idx <- match(g.weight_table$rsid, snp_info$id)
#     snp_pos <- as.integer(snp_info$pos[snp_idx])
#     snp_chrom <- snp_info$chrom[snp_idx]
#     if (any(unique(snp_chrom) != chrom)) {
#         stop(sprintf("'chrom' in weight of %s does not match with the LD reference!", 
#             molecular_id))
#     }
#     wgt.snpinfo <- data.frame(chrom = chrom, id = g.weight_table$rsid, 
#         cm = 0, pos = snp_pos, alt = g.weight_table$eff_allele, 
#         ref = g.weight_table$ref_allele, stringsAsFactors = FALSE)
#     if (verbose) 
#         loginfo("Harmonize weight")
#     harmonized_wgt_res <- harmonize_weights(wgt.matrix, wgt.snpinfo, 
#         snp_info, drop_strand_ambig = drop_strand_ambig)
#     wgt.matrix <- harmonized_wgt_res$wgt.matrix
#     wgt.snpinfo <- harmonized_wgt_res$wgt.snpinfo
#     wgt.matrix <- wgt.matrix[wgt.matrix[, "weight"] != 0, , drop = FALSE]
#     wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]
#     wgt.snp_ids <- intersect(rownames(wgt.matrix), snp_info$id)
#     wgt <- wgt.matrix[match(wgt.snp_ids, rownames(wgt.matrix)), 
#         "weight", drop = FALSE]
#     if (!is.null(top_n_snps)) {
#         wgt <- wgt[order(-abs(wgt[, "weight"])), ]
#         wgt <- head(wgt, top_n_snps)
#         wgt.snp_ids <- rownames(wgt)
#     }
#     if (weight_format == "PredictDB") {
#         if (scale_predictdb_weights) {
#             if (verbose) {
#                 loginfo("Scale weights by variance from LD reference")
#             }
#             wgt_snp_var <- snp_info$variance[match(wgt.snp_ids, 
#                 snp_info$id)]
#             wgt <- wgt * sqrt(wgt_snp_var)
#         }
#     }
#     wgt.snpinfo <- wgt.snpinfo[match(wgt.snp_ids, wgt.snpinfo$id), 
#         ]
#     g.weight_table <- g.weight_table[match(wgt.snp_ids, g.weight_table$rsid), 
#         ]
#     n_wgt <- nrow(wgt)
#     if (n_wgt > 0) {
#         p0 <- min(wgt.snpinfo$pos, na.rm = TRUE)
#         p1 <- max(wgt.snpinfo$pos, na.rm = TRUE)
#         if (!is.null(cov_table)) {
#             if (verbose) {
#                 loginfo("Compute LD of variants in weights (R_wgt) from PredictedDB LD")
#             }
#             g.cov_table <- cov_table[cov_table$GENE == molecular_id, 
#                 ]
#             R_wgt <- get_weight_LD_from_predictdb(g.cov_table, 
#                 g.weight_table, convert_cov_to_cor = TRUE)
#             R_wgt <- R_wgt[wgt.snp_ids, wgt.snp_ids, drop = FALSE]
#         }
#         else {
#             R_wgt <- NULL
#         }
#     }
#     else {
#         p0 <- p1 <- NA
#         wgt <- R_wgt <- NULL
#     }
#     return(list(chrom = chrom, p0 = p0, p1 = p1, wgt = wgt, R_wgt = R_wgt, 
#         molecular_id = molecular_id, weight_name = weight_name, 
#         type = type, context = context, n_wgt = n_wgt))
# }

# harmonize_weights <- function (wgt.matrix,
#                                wgt.snpinfo,
#                                snp_info,
#                                drop_strand_ambig = TRUE) {

#     target_header <- c("chrom", "id", "pos", "alt", "ref")
#     if (!all(target_header %in% colnames(snp_info))) {
#         stop("snp_info needs to contain the following columns: ",
#             paste(target_header, collapse = " "))
#     }
#     wgt.snpinfo <- wgt.snpinfo[match(rownames(wgt.matrix), wgt.snpinfo$id), 
#     ]
#     snpnames <- intersect(wgt.snpinfo$id, snp_info$id)
#     if (length(snpnames) != 0) {
#         snps.idx <- match(snpnames, wgt.snpinfo$id)
#         ld.idx <- match(snpnames, snp_info$id)
#         qc <- allele.qc(wgt.snpinfo[snps.idx, ]$alt, wgt.snpinfo[snps.idx, 
#             ]$ref, snp_info[ld.idx, ]$alt, snp_info[ld.idx, ]$ref)
#         ifflip <- qc[["flip"]]
#         ifremove <- !qc[["keep"]]
#         flip.idx <- snps.idx[ifflip]
#         wgt.snpinfo[flip.idx, c("alt", "ref")] <- wgt.snpinfo[flip.idx, 
#             c("ref", "alt")]
#         wgt.matrix[flip.idx, ] <- -wgt.matrix[flip.idx, ]
#         if (drop_strand_ambig && any(ifremove)) {
#             remove.idx <- snps.idx[ifremove]
#             wgt.snpinfo <- wgt.snpinfo[-remove.idx, , drop = FALSE]
#             wgt.matrix <- wgt.matrix[-remove.idx, , drop = FALSE]
#         }
#     }
#     return(list(wgt.matrix = wgt.matrix, wgt.snpinfo = wgt.snpinfo))
# }

# cols <- c("gene",
#           "rsid",
#           "varID",
#           "ref_allele",
#           "eff_allele",
#           "weight",
#           "pval",
#           "Chromosome")


# get_weight_LD_from_predictdb <- function (g.cov_table, g.weight_table, convert_cov_to_cor = TRUE) {
#     if (convert_cov_to_cor) {
#         g.cov_table <- convert_predictdb_cov_to_cor(g.cov_table)
#     }
#     if (any(grepl("rs", g.cov_table$RSID1))) {
#         g.weight_table$varID <- g.weight_table$rsid
#     }

#     print(sum(!all(g.weight_table$varID %in% unique(c(g.cov_table$RSID1,
#         g.cov_table$RSID2)))))

#     if (!all(g.weight_table$varID %in% unique(c(g.cov_table$RSID1, 
#         g.cov_table$RSID2)))) {
#         stop("Not all variants in weight_table are in cov_table!")
#     }
#     varIDs <- g.weight_table$varID
#     n_wgt <- length(varIDs)
#     if (n_wgt == 0) {
#         return(NULL)
#     }
#     R_wgt <- diag(n_wgt)
#     if (n_wgt > 1) {
#         snp_pairs <- combn(length(varIDs), 2)
#         R_snp_pairs <- apply(snp_pairs, 2, function(x) {
#             match.idx <- which(g.cov_table$RSID1 == varIDs[x[1]] & 
#                 g.cov_table$RSID2 == varIDs[x[2]])
#             if (length(match.idx) == 0) {
#                 match.idx <- which(g.cov_table$RSID1 == varIDs[x[2]] & 
#                   g.cov_table$RSID2 == varIDs[x[1]])
#             }
#             if (length(match.idx) > 0) {
#                 g.cov_table[match.idx, "VALUE"]
#             }
#             else {
#                 NA
#             }
#         })
#         R_wgt[t(snp_pairs)] <- R_snp_pairs
#         R_wgt[t(snp_pairs[c(2, 1), ])] <- R_snp_pairs
#     }
#     rownames(R_wgt) <- g.weight_table$rsid
#     colnames(R_wgt) <- g.weight_table$rsid
#     return(R_wgt)
# }