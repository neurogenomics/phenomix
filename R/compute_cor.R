#' Compute all pairwise trait correlations
#'
#' Computes pairwise correlations for all traits
#'  using either \link[WGCNA]{cor} (faster) or
#' \link[stats]{cor} (slower).
#'
#' @param obj Either a \pkg{Seurat} object or a feature*trait matrix.
#' @param assay If \code{obj} is a \pkg{Seurat} object, which assay to extract
#'  for correlations.
#' @param slot If \code{obj} is a \pkg{Seurat} object, which slot to extract
#'  for correlations. 
#' @param transpose Whether to transpose the matrix first.
#' @param return_obj Whether to return the \pkg{Seurat} object with a new
#' \code{obj@graphs} slot, or to simply return the sparse correlation matrix.
#' @param verbose Print messages.
#' @inheritParams stats::cor
#' @inheritParams WGCNA::cor
#' @inheritParams scKirby::get_obsm
#'
#' @returns A \pkg{Seurat} object with a new \code{obj@graphs} slot,
#' or the sparse correlation matrix.
#'
#' @keywords internal
#' @importFrom Matrix t
#' @importFrom methods as slot
compute_cor <- function(obj,
                        keys = NULL,
                        graph_name = NULL,
                        assay = NULL,
                        slot = NULL,
                        transpose = FALSE,
                        method = "pearson",
                        fill_na = NULL,
                        use = "all.obs",
                        return_obj = TRUE,
                        nThreads = 1,
                        verbose = TRUE) {
    # devoptera::args2vars(compute_cor)
    
    #### Extract relevant matrix ####
    if (!is.null(keys)) {
        mat <- scKirby::get_obsm(
            obj = obj,
            keys = keys,
            verbose = verbose
        )[[1]]
        #### Important! must transpose ####
        transpose <- TRUE
    } else {
        mat <- scKirby::get_x(
            obj = obj,
            assay = assay,
            slot = slot,
            verbose = verbose
        )
        if(is.list(mat)){
            messager(length(mat),"matrices extracted.",
                     "Selecting only the first one:",
                     shQuote(names(mat)[[1]]), v=verbose)
            mat <- mat[[1]]
        }
    }
    if (isTRUE(transpose)) mat <- Matrix::t(mat)
    if (!is.null(fill_na)) mat[is.na(mat)] <- fill_na
    #### Compute corr ####
    if (is_installed(pkg = "WGCNA")) {
        messager("Computing r with WGCNA.", v = verbose)
        cmat <- WGCNA::cor(x = mat,
                           method = method,
                           use = use,
                           nThreads = nThreads, 
                           verbose = verbose)
    } else {
        messager("Computing r with stats.", v = verbose)
        cmat <- stats::cor(x = mat,
                           use = use,
                           method = method)
    }
    #### Convert to sparse graph ####
    cmat <- scKirby::to_graph(cmat)
    if (isTRUE(return_obj) && 
        scKirby::is_class(obj,"seurat")) {
        graph_name <- infer_graph_name(obj = obj, 
                                       graph_name = graph_name,
                                       assay = assay, 
                                       keys = keys, 
                                       ignore_has_graph = TRUE, 
                                       verbose = FALSE)
        messager("Adding new graph to obj:", graph_name, v = verbose)
        obj@graphs[[graph_name]] <- cmat
        return(obj)
    } else {
        messager("Returning sparse correlation matrix.", v = verbose)
        return(cmat)
    }
}
