#' Is package installed
#' 
#' @keywords internal
#' @importFrom utils installed.packages
is_installed <- function(pkg) {
    pkg %in% rownames(utils::installed.packages())
}
