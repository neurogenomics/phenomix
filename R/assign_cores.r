
#' Assign cores 
#' 
#' Automatically infer the number of available cores and
#'  allocate the majority of them to parallel processing.
#'  Reserves a small proportion of cores to avoid slowing or crashing.
#'
#' @export 
#' @importFrom DelayedArray setAutoBPPARAM 
#' @importFrom BiocParallel MulticoreParam 
#' @importFrom parallel detectCores
assign_cores <- function(worker_cores=.90,
                         verbose=T){
  # Enable parallelization of HDF5 functions
  ## Allocate ~10% of your available cores to non-parallelized processes
  worker_cores <- if(is.null(worker_cores)) .90 else worker_cores
  total_cores <- parallel::detectCores()
  if(worker_cores<1){
    reserved_cores <-  ceiling(total_cores*(1-worker_cores))
    workers <- total_cores - reserved_cores
  } else {
    workers <- worker_cores
    reserved_cores <-  total_cores - workers
  }
  message("+ ",workers," core(s) assigned as workers (",reserved_cores," reserved).")
  DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers))
  DelayedArray:::set_verbose_block_processing(verbose)
  return(list(worker_cores=workers,
              reserved_cores=reserved_cores,
              total_cores=total_cores))
}