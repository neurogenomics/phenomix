#' Get a CellTypeDataset
#'
#' Get a CellTypeDataset object.
#' @returns \pkg{CTD} object.
#' 
#' @export
#' @examples 
#' ctd <- get_ctd()
get_ctd <- function(fname="BlueLake2018_FrontalCortexOnly.rds") {
    tmp <- get_data(fname = fname)
    obj <- readRDS(tmp)
    return(obj)
}
