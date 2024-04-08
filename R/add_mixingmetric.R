#' Add mixing metric
#' 
#' Compute the mixing metric for a given reduction and add it to the metadata 
#' of a \link{Seurat} object. See \link[Seurat]{MixingMetric} for more details.
#' @param invert Invert the mixing metric such that higher values 
#' indicate better mixing.
#' @param normalise Normalise the mixing metric to the range [0,1].
#' @inheritParams Seurat::MixingMetric
#' @inheritDotParams Seurat::MixingMetric
#' @export
add_mixingmetric <- function(obj,
                             reduction,
                             new_col=paste0("MixingMetric_",reduction),
                             grouping.var = "source",
                             k = 5,
                             max.k = 300, 
                             dims = NULL,
                             invert = TRUE,
                             normalise = FALSE,
                             ...){
    if(is.null(dims)){
        dims <- seq(ncol(obj[[reduction]]))
    }
    message("Computing mixing metric for reduction: ",reduction,
            " (",max(dims)," dimensions)")
    obj@meta.data[[new_col]] <- Seurat::MixingMetric(
        obj,
        reduction = reduction,
        grouping.var = grouping.var,
        k = k,
        max.k = max.k, 
        dims = dims,
        ...) 
    if(isTRUE(invert)){
        obj@meta.data[[new_col]] <- max.k - obj@meta.data[[new_col]]
    }
    if(isTRUE(normalise)){
        obj@meta.data[[new_col]] <- obj@meta.data[[new_col]] / max.k
    } 
    return(obj)
}