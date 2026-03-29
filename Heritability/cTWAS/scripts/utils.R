################################################################################
################################################################################

prepare_region_info <- function(region_file, CHROM = NULL) {
  region_info <- readRDS(region_file)
  if(!is.null(CHROM)) { region_info <- region_info[region_info$chrom == CHROM, ] }
  return(region_info)
}

################################################################################
################################################################################

str2bool <- function(input_str) {
  if (tolower(input_str) %in% c("true", "t", "1", 1)) {
    return(TRUE)
  } else if (tolower(input_str) %in% c("false", "f", "0", 0)) {
    return(FALSE)
  } else {
    stop("Invalid boolean value: ", input_str)
  }
}

################################################################################
################################################################################

extract_base_dir_label <- function(path) {
  child_dir <- basename(path)
  return(gsub(".vcf|.gz|.csv|.tsv|.rds|.RDS", "", child_dir))
}

################################################################################
################################################################################

save_RDS <- function(data, prefix) {
  
  dir.create(dirname(prefix),
             recursive = TRUE, showWarnings = FALSE)

  saveRDS(data, glue("{prefix}.rds"))

}

################################################################################
################################################################################

paste_outpath <- function(path, folder_name, add_trailing_slash = TRUE) {
  path_clean <- sub("[/\\\\]+$", "", path)
  outpath <- file.path(path_clean, folder_name)
  outpath <- if (add_trailing_slash) paste0(outpath, .Platform$file.sep) else outpath
  return(outpath)
}

################################################################################
################################################################################

generate_outpath <- function(path, folder_name) {
  outpath <- paste_outpath(path, folder_name)
  dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
  return(outpath)
}

################################################################################
################################################################################

is_scalar <- function(x) is.atomic(x) && length(x) == 1

################################################################################
################################################################################

path_parser <- function(paths, key = NULL) {
  # Case 1: YAML file
  if (is_scalar(paths) && tools::file_ext(paths) %in% c("yaml", "yml")) {
    paths <- yaml::read_yaml(paths)

    if (!is.null(key)) paths <- paths[[key]]
    if (is.character(paths)) paths <- as.list(paths)
    if (!is.list(paths)) stop("YAML must contain a list of file paths")

    # Case 2: anything else that's not already a list
  } else if (is_scalar(paths) && !is.list(paths)) {
    paths <- list(paths)

    # Optional type check for lists
  } else if (!all(sapply(paths, is.character))) {
    stop("All elements must be strings or paths")
  }

  return(paths)
}

################################################################################
################################################################################