#' Run LDA
#'
#' Run Linear Discriminant Analysis (LDA).
#' \emph{WARNING:} Takes a very long time to run on medium to large datasets.
#'
#'
#' Uses \link[stats]{prcomp}.
#' @source \url{https://towardsdatascience.com/linear-discriminant-analysis-lda-101-using-r-6a97217a55a6}
#'
#' @param mat Matrix to run PCA on.
#' @param transpose Whether to transpose the matrix first.
#' @param ... Additional parameters passed to \link[MASS]{lda}.
#' @inheritParams MASS::lda
#'
#' @importFrom Matrix t
#' @importFrom MASS lda
#'
#' @examples
#' mat <- as.matrix(iris[, seq(1, 4)])
#' grouping <- iris$Species
#' @export
run_lda <- function(mat,
                    grouping,
                    transpose = TRUE,
                    center = TRUE,
                    scale. = FALSE,
                    rank. = NULL,
                    ...) {
    if (transpose) {
        mat <- Matrix::t(mat)
    }
    lda_res <- MASS::lda(
        x = mat,
        grouping = grouping
    )
    return(lda_res)
}
