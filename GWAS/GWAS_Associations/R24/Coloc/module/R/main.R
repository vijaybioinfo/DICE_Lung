METHODS.ROOT <- list("abf" = "abf/main.R")

add.prefix.to.methods.root <- function(prefix){
    NEW.ROOTS <- paste0(prefix, "/", METHODS.ROOT)
    names(NEW.ROOTS) <- names(METHODS.ROOT)
    assign("METHODS.ROOT", as.list(NEW.ROOTS), envir = .GlobalEnv)
}


load.system <- function(root){
    ######## Remember Working Directory ############
    CWD <- getwd()
    ######## Source System ##########
    source(root, chdir=TRUE)
    ######## Back to original working directory ###########
    setwd(CWD)
}


run.interface <- function(method, params){
    ########### LOAD METHOD #############
    load.system(METHODS.ROOT[[method]])

    ########### GET METHOD PARAMS ############
    #args <- set.params(params, args.default.list)
    
    ########### RUN METHOD ############
    run.from.cli(params)
}


#####################################################
####################### CLI #########################
#####################################################
if (sys.nframe() == 0L){
    options <- commandArgs(trailingOnly = FALSE)
    script <- sub("--file=", "", options[grep("--file=", options)])
    add.prefix.to.methods.root(dirname(script))

    args <- commandArgs(trailingOnly = TRUE)

    if (!(args[1] %in% names(METHODS.ROOT))){
        stop(paste("Please select one of the following programs {", names(METHODS.ROOT), "}"))
    }

    run.interface(args[1], args[-1])
}
