#' Get a CellTypeDataset
#'
#' Get a CellTypeDataset object.
#' @inheritParams piggyback::pb_download
#' @returns \pkg{CTD} object.
#' 
#' @export
#' @examples 
#' ctd <- get_ctd()
get_ctd <- function(file="BlueLake2018_FrontalCortexOnly.rds") {
    tmp <- get_data(file = file)
    obj <- readRDS(tmp)
    return(obj)
}
