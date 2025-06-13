.here <- function(){
    options <- commandArgs(trailingOnly = FALSE)
    script <- sub("--file=", "", options[grep("--file=", options)])
    return(dirname(script))
}

.source.external <- function(path){
    source(paste0(.here(), "/", path), chdir=TRUE)
}

if (sys.nframe() == 0L){
    if (exists("snakemake")){
        source("module/R/main.R", chdir=TRUE)
        run.interface(snakemake@params$method, snakeamke@params)

    } else {
        .source.external("../../module/R/main.R")
        add.prefix.to.methods.root(paste0(.here(), "/", "../../module/R"))

        args <- commandArgs(trailingOnly = TRUE)
        run.interface(args[1], args[-1])
    }
}