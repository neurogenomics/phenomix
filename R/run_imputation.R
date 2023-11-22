#' Run imputation
#' 
#' Run matrix imputation with various methods.
#' @param X Matrix to impute.
#' @param Xnet Network coefficients to use for imputation.
#' @inheritParams ADImpute::Impute
#' @inheritDotParams ADImpute::Impute
#' @return Imputed matrix.
#' 
#' @export
#' @examples
#' X <- scKirby::example_obj("matrix")
#' Xnet <- ADImpute::network.coefficients
#' imputed <- run_imputation(X, Xnet)
run_imputation <- function(X,
                           Xnet=ADImpute::network.coefficients,
                           do = c("Network"), 
                           cores = 1,
                           ...){
    imputed <- ADImpute::Impute(data = X, 
                                do = do, 
                                cores = cores,
                                net.coef = Xnet, 
                                ...
                      )
    
    return(imputed)
}