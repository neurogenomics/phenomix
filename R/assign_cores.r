#' Assign cores
#'
#' Automatically infer the number of available cores and
#'  allocate the majority of them to parallel processing.
#'  Reserves a small proportion of cores to avoid slowing or crashing.
#' @inheritParams BiocParallel::MulticoreParam
#'
#' @keywords internal
#' @importFrom DelayedArray setAutoBPPARAM
#' @importFrom BiocParallel MulticoreParam
#' @importFrom parallel detectCores
assign_cores <- function(workers = .90,
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
    messager("+",workers,"core(s) assigned as workers",
             paste0("(",reserved_cores, " reserved)."),v=verbose)
    BPPARAM <- BiocParallel::MulticoreParam(workers = workers,
                                            progressbar = progressbar)
    DelayedArray::setAutoBPPARAM(BPPARAM = BPPARAM)
    # DelayedArray:::set_verbose_block_processing(verbose)
    return(BPPARAM)
}
