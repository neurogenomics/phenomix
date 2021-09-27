#' Create a \code{DimReducObject} from DeGAs
#'
#' Create a \pkg{Seurat} \code{DimReducObject}
#' from DeGAs embeddings and loadings.
#'
#' @param reconstruction Output of \code{reconstruct_matrix}
#'  with embeddings and loadings.
#' @param assay Assay to use.
#' @param key Key to use (name of embedding).
#'
#' @keywords internal
#' @importFrom Seurat CreateDimReducObject
create_DimReducObject <- function(reconstruction,
                                  assay,
                                  key) {
    DR <- Seurat::CreateDimReducObject(
        embeddings = reconstruction$embeddings,
        loadings = reconstruction$loadings,
        stdev = if (is.null(reconstruction$stdev)) {
            numeric()
        } else {
            as.numeric(reconstruction$stdev)
        },
        assay = assay,
        key = key
    )
    return(DR)
}
