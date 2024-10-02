#' Run correlation
#' 
#' Run pairwise correlations on a matrix.
#' @export
#' @param X Matrix to compute correlations
#' @param t Transpose the matrix.
run_cor <- function(X,
                    t=FALSE){
    if(t) X <- Matrix::t(X)
    message("Computing correlations.")
    ## Computing dist first improve dynamic range of similarity
    ## Much faster than stats::dist 
    # Xd <- Rfast::Dist(as.matrix(X))
    # rownames(Xd) <- colnames(X)
    # colnames(Xd) <- colnames(X)
    # Xd
    ## Convert to similarities (same approach as proxy::simil)
    # Xd <- 1/(1 + Xd)
    # proxy::simil(d)|>as.matrix() 
    WGCNA::cor(X, cosine = TRUE)
}