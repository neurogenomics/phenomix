#' Run PCA
#'
#' Run Principal Components Analysis (PCA).
#'
#' Uses \link[stats]{prcomp}.
#'
#' @param mat Matrix to run PCA on.
#' @param transpose Whether to transpose the matrix first.
#' @param ... Additional parameters passed to \link[stats]{prcomp}.
#' @inheritParams stats::prcomp
#'
#' @importFrom Matrix t
#' @importFrom stats prcomp
#'
#' @export
run_pca <- function(mat,
                    transpose = TRUE,
                    center = TRUE,
                    scale. = FALSE,
                    rank. = NULL,
                    ...) {
    if (transpose) {
        mat <- Matrix::t(mat)
    }
    pca <- stats::prcomp(
        x = mat,
        center = center,
        scale. = scale.,
        rank. = rank.,
        ...
    )
    return(pca)
}
