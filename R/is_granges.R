#' Is an object of class GRanges
#'
#' Check whether an object is of the class \link[GenomicRanges]{GRanges}.
#' @param obj Any R object.
#' @returns Boolean.
#' @export
#' @importFrom methods is
#' @examples 
#' bool <- is_granges(mtcars)
is_granges <- function(obj) {
    methods::is(obj, "GRanges")
}
