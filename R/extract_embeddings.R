#' Extract embeddings
#'
#' Extract embeddings from various object types.
#'
#' @param obj Object of class:
#' \itemize{
#' \item{\code{DimReduc}}
#' \item{\code{Seurat}}
#' \item{\code{Seurat}}
#' }
#' @param reduction Which reduction to extract from.
#' Only used if \code{obj} is of class \code{Seurat}.
#' @param verbose Print messages.
#'
#' @export
#' @importFrom methods is
#' @importFrom SeuratObject Reductions
extract_embeddings <- function(obj,
                               reduction = NULL,
                               verbose = TRUE) {

    # obj <- DEGAS_contributionGene; obj <- scNLP::pseudo_seurat
    if (methods::is(obj, "DimReduc")) {
        messager("Extracting embeddings from DimReduc.", v = verbose)
        DR <- obj
        embeddings <- DR@cell.embeddings
    } else if (methods::is(obj, "Seurat")) {
        if (is.null(reduction)) reduction <- Seurat::Reductions(obj)[1]
        messager("Extracting embeddings from Seurat reduction:",
            reduction,
            v = verbose
        )
        DR <- obj@reductions[[reduction]]
        embeddings <- DR@cell.embeddings
    } else if (methods::is(obj, "prcomp")) {
        messager("Extracting embeddings from prcomp.", v = verbose)
        embeddings <- obj$x
    } else if (methods::is(obj, "list") &
        methods::is(obj, "AssayData")) {
        messager("Extracting embeddings from list", v = verbose)
        if ("embedding" %in% names(obj)) {
            embeddings <- obj$embedding
        } else if ("u" %in% names(obj)) {
            embeddings <- obj$u
        } else if ("U" %in% names(obj)) {
            embeddings <- obj$U
        } else {
            stop("No embeddings could be identified.")
        }
    } else {
        stop("No embeddings could be identified.")
    }
    return(embeddings)
}
