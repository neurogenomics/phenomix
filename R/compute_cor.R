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
#' @param reduction If not \code{NULL}, extracts the embedding from the
#' specified reduction (e.g. "pca") and uses that to
#' compute trait-trait similarities instead.
#' @param transpose Whether to transpose the matrix first.
#' @param return_obj Whether to return the \pkg{Seurat} object with a new
#' \code{obj@graphs} slot, or to simply return the sparse correlation matrix.
#' @inheritParams stats::cor
#'
#' @returns A \pkg{Seurat} object with a new \code{obj@graphs} slot,
#' or the sparse correlation matrix.
#'
#' @export
#' @importFrom Matrix t
#' @importFrom methods as
compute_cor <- function(obj,
                        assay = NULL,
                        slot = NULL,
                        reduction = NULL,
                        transpose = FALSE,
                        method = "pearson",
                        return_obj = TRUE,
                        verbose = TRUE) {
    #### Extract relevant matrix ####
    if (!is.null(reduction)) {
        mat <- extract_embeddings(
            obj = obj,
            reduction = reduction,
            verbose = verbose
        )
    } else {
        mat <- extract_matrix(
            obj = obj,
            assay = assay,
            slot = slot,
            verbose = verbose
        )
    }
    if (transpose) mat <- Matrix::t(mat)
    #### Compute corr ####
    if (is_installed(pkg = "WGCNA")) {
        messager("Computing r with WGCNA.", v = verbose)
        cmat <- WGCNA::cor(mat, method = method)
    } else {
        messager("Computing r with stats", v = verbose)
        cmat <- stats::cor(mat, method = method)
    }
    cmat <- methods::as(methods::as(cmat, "sparseMatrix"), "Graph")
    if (return_obj & is_seurat(obj)) {
        graph_name <- if (is.null(assay)) "cor" else paste0(assay, "_", "cor")
        messager("Adding new graph to obj:", graph_name, v = verbose)
        obj@graphs[[graph_name]] <- cmat
        return(obj)
    } else {
        messager("Returning sparse correlation matrix.", v = verbose)
        return(cmat)
    }
}
