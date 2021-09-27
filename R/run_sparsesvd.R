#' Run sparse SVD
#'
#' Run sparse Singular Value Decomposition (sparseSVD).
#'
#' Uses \link[sparsesvd]{sparsesvd}.
#'
#' @param mat Matrix to run sparseSVD on.
#' @param transpose Whether to transpose the matrix first.
#' @param add_names Add colnames and rownames to embeddings and loadings.
#' @inheritParams sparsesvd::sparsesvd
#'
#' @importFrom Matrix t
#' @importFrom sparsesvd sparsesvd
#'
#' @export
run_sparsesvd <- function(mat,
                          transpose = TRUE,
                          add_names = TRUE,
                          rank = 0L,
                          tol = 1e-15,
                          kappa = 1e-6) {
    if (transpose) {
        mat <- Matrix::t(mat)
    }
    mat <- Matrix::Matrix(mat, sparse = TRUE)
    ssvd <- sparsesvd::sparsesvd(
        M = mat,
        rank = rank,
        tol = tol,
        kappa = kappa
    )
    if (add_names) {
        ssvd_names <- paste0("SSVD.", seq(1, ncol(ssvd$u)))
        ssvd$u <- ssvd$u %>%
            `colnames<-`(ssvd_names) %>%
            `row.names<-`(rownames(mat))
        ssvd$v <- ssvd$v %>%
            `colnames<-`(ssvd_names) %>%
            `row.names<-`(colnames(mat))
        ssvd$d <- setNames(ssvd$d, ssvd_names)
    }
    return(ssvd)
}
