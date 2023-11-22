#' Set cores
#'
#' Assign cores automatically for parallel processing, while reserving some.
#' 
#' @param workers Number (>1) or proportion (<1) of worker cores to use.
#' @param verbose Print messages.
#' @param progressbar logical(1) Enable progress bar 
#' (based on \code{plyr:::progress_text}).
#'  Enabling the progress bar changes the default value of tasks to 
#'  \code{.Machine$integer.max}, so that progress is reported for 
#'  each element of X.
#' @return List of core allocations.
#'
#' @importFrom DelayedArray setAutoBPPARAM
#' @importFrom parallel detectCores
#' @keywords internal
set_cores <- function(workers = .90,
                         progressbar = TRUE,
                         verbose = TRUE) {
    # Enable parallelization of HDF5 functions
    ## Allocate ~10% of your available cores to non-parallelized processes
    workers <- if (is.null(workers)) .90 else workers
    total_cores <- parallel::detectCores()
    if (workers < 1) {
        reserved_cores <- ceiling(total_cores * (1 - workers))
        workers <- total_cores - reserved_cores
    } else {
        workers <- workers
        reserved_cores <- total_cores - workers
    }
    messager(workers, "core(s) assigned as workers",
             paste0("(",reserved_cores, " reserved)."),
             v = verbose
    )
    #### Handle Windows ####
    if (.Platform$OS.type == "windows") {
        params <- BiocParallel::SnowParam(workers = workers,
                                          progressbar = progressbar)
    } else {
        params <- BiocParallel::MulticoreParam(workers = workers, 
                                               progressbar = progressbar)
    }
    DelayedArray::setAutoBPPARAM(params)
    #### Not allowed to use internal functions ####
    # DelayedArray:::set_verbose_block_processing(verbose)
    return(list(
        workers = workers,
        reserved_cores = reserved_cores,
        total_cores = total_cores,
        params = params
    ))
}
