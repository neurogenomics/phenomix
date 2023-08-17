#' Is package installed
#' 
#' @keywords internal
#' @importFrom utils installed.packages
is_installed <- function(pkg) {
    ii <- pkg %in% rownames(utils::installed.packages())
    if(isFALSE(ii)){
        messager("Warning: Must install",shQuote(pkg),"to use this feature.")
    }
    return(ii)
}
