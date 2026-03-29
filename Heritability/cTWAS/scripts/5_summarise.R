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

paste_path <- function(name, path) { glue("{path}{name}") }

# Function to make compatible with snakemake but not
# break original modality of providing dir name
# rather than list of paths
extract_path <- function(files_path) {

  # If single directory → expand
  if (length(files_path) == 1 && dir.exists(files_path)) {

    files <- list.files(
      files_path,
      pattern = "\\.RDS$",
      full.names = TRUE
    )
    print("One file passed")

  } else {

    # Assume vector of file paths
    files <- unlist(files_path)
    print("Multiple files passed")
    print(class(files))
    # Basic validation

    if (!all(file.exists(files))) {
      missing <- files[!file.exists(files)]
      stop("Missing files:\n", paste(missing, collapse = "\n"))
    }

  }

  # Keep only RDS files
  files <- files[grepl("\\.RDS$", files)]

  if (length(files) == 0) {
    stop("No .RDS files found.")
  }

  names(files) <- tools::file_path_sans_ext(basename(files))
  print(names(files))

  return(files)
}

# extract_path <- function(files_path) {

  
#   files_path <- list.files(files_path,
#                            pattern = "\\.RDS$",
#                            full.names = TRUE)

#   names(files_path) <- gsub("\\.RDS",
#                             "",
#                             basename(files_path))

#   return(files_path)

# }

################################################################################
################################################################################

extract_gwas_n <- function(metatable_path, study) {
  table <- read.csv(metatable_path, sep=",")
  gwas_n <- table[table$Sample.name == study, ]$sample_size
  return(as.integer(gwas_n))
}

################################################################################
################################################################################

prepare_gene_annotation <- function(ens_db, finemap_res) {
  finemap_gene_res <- subset(finemap_res, group != "SNP")
  gene_ids <- unique(finemap_gene_res$molecular_id)
  gene_annot <- ctwas::get_gene_annot_from_ens_db(ens_db, gene_ids)
  return(gene_annot)
}

################################################################################
################################################################################

add_gene_annotation <- function(finemap_res, snp_map, gene_annotation) {

  colnames(gene_annotation)[colnames(gene_annotation) == "gene_id"] <- "molecular_id"

  finemap_res <- ctwas::anno_finemap_res(finemap_res,
                                         snp_map = snp_map,
                                         mapping_table = gene_annotation,
                                         add_gene_annot = TRUE,
                                         map_by = "molecular_id",
                                         drop_unmapped = FALSE,
                                         add_position = TRUE,
                                         use_gene_pos = "mid")

  return(finemap_res)
}

extract_finemap_df <- function(name, paths) {
  df <- readRDS(paths[[name]])
  df$pval <- ctwas::z2p(df$z)
  df$greater_group <- gsub(".RDS|.rds", "", name)
  return(df)
}


prepare_finemap_res <- function(finemap_path, snp_map) {

  paths <- extract_path(finemap_path)

  finemap_res_list <- lapply(names(paths),
                             extract_finemap_df,
                             paths = paths)

  # finemap_res <- dplyr::bind_rows(finemap_res_list)
  finemap_res <- data.table::rbindlist(finemap_res_list, use.names = TRUE, fill = TRUE)

  finemap_res <- data.frame(finemap_res)
  
  if (dim(finemap_res)[1] == 0) {
    message("Empty Dataframe finemap_res dataframe for ", finemap_path)
    return(finemap_res)
  }

  gtf <- prepare_gene_annotation(EnsDb.Hsapiens.v86, finemap_res)
  finemap_res <- add_gene_annotation(finemap_res, snp_map, gtf)

  return(finemap_res)
}

################################################################################
################################################################################

extract_heritability <- function(ctwas_parameters, group_label) {
  data <- data.table::data.table(
    category = names(ctwas_parameters$prop_heritability),
    percentage = ctwas_parameters$prop_heritability,
    group_label = group_label
  )
  data$percentage_label <- paste0(round(data$percentage * 100), "%")
  return(data)
}

internal_lapply_herit_fun <- function(group_label, obj, gwas_n) {
  param <- readRDS(obj[[group_label]])

  ctwas_parameters <- ctwas::summarize_param(param,
                                             gwas_n,
                                             enrichment_test = "fisher")

  data <- extract_heritability(ctwas_parameters, gsub("\\.RDS", "", group_label))

  return(data)
}

herit_table <- function(input_path, gwas_n) {

  paths <- extract_path(input_path)

  dfs <- lapply(names(paths),
                FUN = internal_lapply_herit_fun,
                obj = paths,
                gwas_n = gwas_n)

  

  # df <- dplyr::bind_rows(df)
  df <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)

  return(df)
}

################################################################################
################################################################################

generate_frame <- function(name, param_sub_group_obj) {
  vect <- param_sub_group_obj[[name]]

  df <- data.table::data.table(
    name     = names(vect),
    score    = unname(vect),
    category = name
  )

  # # Ensure name column exists even if names(v) is NULL
  # if (is.null(df$name)) df[, name := NA_character_]

  return(df)
}

summarise_parameters_single <- function(group_label, obj, gwas_n) {
  print(group_label)
  param <- readRDS(obj[[group_label]])
  if ("empty_list_flag" %in% unlist(param) | is.null(param)) {
    return(data.table(
      "name" = character(),
      "score" = numeric(),
      "category" = character(),
      "group_level" = character(),
      "label" = character()
    ))
  }

  ctwas_parameters <- ctwas::summarize_param(
    param,
    gwas_n,
    enrichment_test = "fisher"
  )

  df <- data.table::rbindlist(
    lapply(names(ctwas_parameters),
           FUN = generate_frame,
           param_sub_group_obj = ctwas_parameters),
    use.names = TRUE,
    fill = TRUE
  )

  # gene names present in df (excluding SNP and NA)
  gene_names <- df[!is.na(name) & name != "SNP", unique(name)]

  # gene names present in df (excluding SNP and NA)
  gene_names <- df[!is.na(name) & name != "SNP", unique(name)]

  # if (length(gene_names) > 1L) {
  #   stop(
  #     sprintf(
  #       "Expected at most one tissue-context (non-SNP name) in df, found %d: %s",
  #       length(gene_names),
  #       paste(gene_names, collapse = ", ")
  #     )
  #   )
  # }
  
  # Add explicit totals for SNP and gene (if a gene name exists)
  tot_base <- df[category == "total_pve" & !is.na(score)]
  tot_snp <- data.table::copy(tot_base)
  tot_snp[, name := "SNP"]
  tot_gene <- NULL

  # total_pve rows expanded for each gene name
  df <- data.table::rbindlist(
    list(
      df[!is.na(name)],
      tot_snp,
      tot_gene
    ),
    use.names = TRUE,
    fill = TRUE
  )

  df[, group_label := gsub("\\.RDS$", "", group_label)]

  df
  return(df)
}

# summarise_parameters_iterative <- function(input_path, gwas_n) {


#   paths <- extract_path(input_path)
#   print(paths)

#   dfs <- lapply(names(paths),
#                 FUN = summarise_parameters_single,
#                 obj = paths,
#                 gwas_n = gwas_n)

#   # df <- dplyr::bind_rows(dfs)
#   dfl <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)

#   dfl[, label :=
#     fcase(
#       name == "SNP", paste(name, group_label, sep = "_"),
#       default = paste0("gene_", group_label)
#     )
#   ]


#   dfw <- dcast(dfl, label ~ category, value.var = "score")

#   return(list(`long` = df, `wide` = dfw))
# }

summarise_parameters_iterative <- function(param_path, gwas_n) {

  paths <- extract_path(param_path)

  if (length(paths) == 0) stop("No input RDS files found in: ", param_path)

  dfs <- lapply(
    names(paths),
    FUN = summarise_parameters_single,
    obj = paths,
    gwas_n = gwas_n
  )

  dfl <- data.table::rbindlist(dfs, use.names = TRUE, fill = TRUE)

  if (dim(dfl)[1] == 0 ) { 
    message("Empty params data.table for ", paths)
    return (
    list(long = dfl, wide = dfl)) 
    }

  dfl[, label := data.table::fcase(
    name == "SNP", paste(name, group_label, sep = "_"),
    default = paste0("gene_", group_label)
  )]

  dfw <- data.table::dcast(dfl, label ~ category, value.var = "score")

  return( list(long = dfl, wide = dfw) )
}

################################################################################
################################################################################

# combined_pip_by_type <- combine_gene_pips(ctwas_obj[["no_model_eqtl_Adipose"]][["finemap_res"]],
#                                           # mapping_table = gtf,
#                                           group_by = "molecular_id",
#                                           by = "type",
#                                           method = "combine_cs",
#                                           filter_cs = FALSE,
#                                           include_cs_id = TRUE)

################################################################################
################################################################################
save_files <- function(List, prefix) {

  # Create output directory if it does not exist
  dir.create(prefix, recursive = TRUE, showWarnings = FALSE)

  # Save each element as TSV
  lapply(names(List), function(x) {

    file_path <- file.path(prefix, paste0(x, ".tsv"))
  
    print(head(file_path))
    print(List[[x]])
  
    write.table(
      List[[x]],
      file = file_path,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  })
}
################################################################################
################################################################################

# find_gwas_n <- function(metatable_path, study, gwas_n) {

#   if (!is.null(metatable_path) && !is.null(study)) {
#     gwas_n <- extract_gwas_n(metatable_path, study)

#   } else if (!is.null(gwas_n)) {
#     gwas_n <- gwas_n

#   } else {
#     stop(
#       "WARNING - Must provide EITHER metatable_path AND study name OR gwas_n"
#     )
#   }
#   return(gwas_n)
# }

################################################################################
################################################################################

# main <- function(finemap_path,
#                  param_path,
#                  screen_res_path,
#                  susie_alpha_res_path,
#                  z_gene_path,
#                  snp_map_path,
#                  gwas_n)

main <- function(finemap_path,
                 param_path,
                 snp_map_path,
                 gwas_n) {

  snp_map <- readRDS(snp_map_path)

  finemap_res <- prepare_finemap_res(finemap_path, snp_map)
  print(head(finemap_res))
  print(class(finemap_res))

  parameter_df_list <- summarise_parameters_iterative(param_path, gwas_n)

  heritability_table <- herit_table(param_path, gwas_n)

  # finemap_res_ss <- subset(finemap_res, group != "SNP" & susie_pip > 0.8 & !is.na(cs))
  message("Files processed")
  return( list(`finemap_res` = finemap_res,
               `finemap_res_filt_0.7` = finemap_res[finemap_res$susie_pip > 0.7, ],
               `parameters_long` = parameter_df_list[['long']],
               `parameters_wide` = parameter_df_list[['wide']],
               `heritability_table` = heritability_table)
  )
}

################################################################################
################################################################################

if (sys.nframe() == 0) {

  script_dir <- get_script_dir()
  source(file.path(script_dir, "utils.R"))

  suppressPackageStartupMessages({
    library("glue")
    library("dplyr")
    library("ctwas")
    library("enrichR")
    library("ggplot2")
    library("argparse")
    library("data.table")
    library("EnsDb.Hsapiens.v86")
  })

  ################################################################################
  ################################################################################

  parser <- ArgumentParser(description = "Prepare reference files for cTWAS")

  parser$add_argument("--finemap_path",
                      nargs="+",
                      help = "",
                      required = TRUE)

  parser$add_argument("--param_path",
                      nargs="+",
                      help = "",
                      required = TRUE)

  parser$add_argument("--snp_map_path",
                      help = "",
                      required = FALSE)

  parser$add_argument("--gwas_n",
                      help = "",
                      required = FALSE,
                      default = NULL,
                      type = "integer")

  parser$add_argument("--outpath",
                      help = "outpath to .RDS files",
                      required = TRUE)

  ################################################################################
  ################################################################################

  args <- parser$parse_args()

  # gwas_n <- find_gwas_n( metatable_path = args$metatable_path,
  #                        study = args$study,
  #                        gwas_n = args$gwas_n )

  output <- main( finemap_path = args$finemap_path,
                  param_path = args$param_path,
                  snp_map_path = args$snp_map_path,
                  gwas_n = args$gwas_n )

  save_files(output, args$outpath)

  ################################################################################
  ################################################################################
}

  # output <- main( finemap_path = args$finemap_path,
  #                 param_path = args$param_path,
  #                 screen_res_path = args$screen_res_path,
  #                 susie_alpha_res_path = args$susie_alpha_res_path,
  #                 z_gene_path = args$z_gene_path,
  #                 snp_map_path = args$snp_map_path,
  #                 study = args$study,
  #                 gwas_n = args$gwas_n )

  # parser$add_argument("--screen_res_path",
  #                     help = "",
  #                     required = TRUE)

  # parser$add_argument("--susie_alpha_res_path",
  #                     help = "",
  #                     required = TRUE)

  # parser$add_argument("--z_gene_path",
  #                     help = "",
  #                     required = TRUE)
