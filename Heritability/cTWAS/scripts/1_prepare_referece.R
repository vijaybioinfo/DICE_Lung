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

prepare_region_metatable_UK_biobank <- function(genome_version,
                                                region_info,
                                                ld_dir) {

  print(
    paste0("Processing for chromosmoe - ",
           unique(region_info$chrom))
  )

  ld_filestem <- sprintf("ukb_%s_0.1_chr%s.R_snp.%s_%s",
                         genome_version,
                         region_info$chrom,
                         region_info$start,
                         region_info$stop)

  region_metatable <- region_info
  region_metatable$LD_file <- file.path(ld_dir, paste0(ld_filestem, ".RDS"))
  region_metatable$SNP_file <- file.path(ld_dir, paste0(ld_filestem, ".Rvar"))
  region_metatable$label <- sub("\\.\\d+_\\d+\\.RDS$",
                                "",
                                basename(region_metatable$LD_file))
  return(region_metatable)
}

################################################################################
################################################################################

extract_available_ld_dir_paths_ <- function(genome_version,
                                            region_info,
                                            ld_dir) {

  regex_pattern <- ("^(?i)(chr)?(1[0-9]|2[0-2]|[1-9]|X|Y)$")
  sub_dirs <- list.files(ld_dir)
  sub_dirs <- sub_dirs[grep(regex_pattern, sub_dirs)]
  ld_dir_paths <- paste_outpath(ld_dir, sub_dirs, add_trailing_slash = FALSE)
  return(ld_dir_paths)

}

################################################################################
################################################################################

prepare_res <- function(region_info, genome_version, ld_dir) {

  region_metatable <- prepare_region_metatable_UK_biobank(genome_version,
                                                          region_info,
                                                          ld_dir)

  res <- create_snp_LD_map(region_metatable)
  res$label <- unique(region_metatable$label)

  invisible(res)
}

################################################################################
################################################################################

res_signle_chr <- function(sargs) {
  region_info <- prepare_region_info(sargs$region_path, sargs$chrom)
  region_metatable <- prepare_region_metatable_UK_biobank(sargs$genome_version,
                                                          region_info,
                                                          sargs$ld_dir)

  res <- create_snp_LD_map(region_metatable)
  res$label <- unique(region_metatable$label)
  return(res)
}

res_multi_chr <- function(margs) {
  region_info <- prepare_region_info(margs$region_path)

  ld_dir_paths <- extract_available_ld_dir_paths_(
                                    margs$genome_version,
                                    margs$region_info,
                                    margs$ld_dir)

  res_dfs <- parallel::mclapply(ld_dir_paths, function(path) {

    splt <- unlist(strsplit(path, .Platform$file.sep))
    chr <- splt[[length(splt)]]

    region_info <- region_info[region_info$chrom == chr, ]

    df <- prepare_region_metatable_UK_biobank(
                        margs$genome_version,
                        region_info,
                        path)
                                              # chromosome = chr)
    data.frame(df)}, mc.cores = 12)

  region_metatable <- do.call(rbind, res_dfs)

  res <- create_snp_LD_map(region_metatable)

  res$label <- sub("_chr.*$", "_chrom_all", region_metatable$label[[1]])
  return(res)
}

################################################################################
################################################################################

save_res <- function(res, outpath, label) {
  sep <- .Platform$file.sep
  decorator <- function(sd) {dir.create(glue("{outpath}{sep}{sd}{sep}"),
                                        recursive=TRUE, showWarnings=FALSE)}

  subdirs <- c("region_info", "snp_map", "LD_map")
  lapply(subdirs, decorator)

  saveRDS(res$region_info, glue("{outpath}{sep}region_info{sep}{label}.rds"))
  saveRDS(res$snp_map, glue("{outpath}{sep}snp_map{sep}{label}.rds"))
  saveRDS(res$LD_map, glue("{outpath}{sep}LD_map{sep}{label}.rds"))
}

################################################################################
################################################################################
################################################################################
################################################################################

main <- function(region_path, ld_dir, genome_version, chrom, outpath) {

  sargs <- list(`region_path` = region_path,
                `ld_dir` = ld_dir,
                `genome_version` = genome_version,
                `chrom` = chrom)

  margs <- list(`region_path` = region_path,
                `ld_dir` = ld_dir,
                `genome_version` = genome_version)

  res <- if (!is.null(chrom)) res_signle_chr(sargs) else res_multi_chr(margs)
  ######################################################
  save_res(res, args$outpath, res$label)

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

  parser$add_argument("--region_path", help = "region path", required = TRUE)

  parser$add_argument("--ld_dir", help = "LD dir", required = TRUE)

  parser$add_argument("--genome_version", help = "genome_version",
                      required = FALSE,
                      default="b38")

  parser$add_argument("--chromosome", type = "integer",
                      help = "chromosme, for subsetting", required = FALSE,
                      default = NULL)

  parser$add_argument("--outpath", help = "outpath to .RDS files",
                      required = TRUE)

  ######################################################

  args <- parser$parse_args()

  ######################################################

  main(region_path = args$region_path,
       ld_dir = args$ld_dir,
       genome_version = args$genome_version,
       chrom = args$chromosome,
       outpath = args$outpath)
}