source.external <- function(path){
    ################ get CWD ####################
    CWD <- getwd()

    ################ Sourcing ###################
    options <- commandArgs(trailingOnly = FALSE)
    script <- sub("--file=", "", options[grep("--file=", options)])
    directory <- dirname(script)
    other.name <- file.path(directory, path)

    source(other.name, chdir=T)

    ################ Back to CWD ################
    setwd(CWD)
}


if (sys.nframe() == 0L){
    if (exists("snakemake")){
        source("module/R/main.R", chdir=T)
        add.prefix.to.methods.root("module/R")
        run.interface(snakemake@params[[1]], snakemake@params[[2]])
    }
}