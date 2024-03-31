#' Add cluster colors
#' 
#' Add or extract cluster colors from a \link{Seurat} object.
#' @export
#' @examples
#' obj <- get_HPO()
#' obj2 <- add_cluster_colors(obj)
add_cluster_colors <- function(obj, 
                               cluster_col = "seurat_clusters",
                               color_col = "cluster_colors",
                               preferred_palettes = "kovesi.cyclic_mygbm_30_95_c78",
                               force_new=FALSE,
                               return_dict=FALSE){
    if(color_col %in% colnames(obj@meta.data) && 
       isFALSE(force_new)){
        d <- obj@meta.data[,c(cluster_col,color_col)]|>unique()|>
            data.table::setorderv(cluster_col)
        messager("Cluster colors already exist in obj:",color_col)
        if(isTRUE(return_dict)) {
            return(
                stats::setNames(d[[color_col]],
                                d[[cluster_col]])
            )
        } else {
            return(obj)
        }
    } else {
        cluster_colors <- KGExplorer::map_colors(
            dat = obj@meta.data,
            columns = cluster_col,
            as = "dict",
            preferred_palettes = preferred_palettes)[[1]]
        if(isTRUE(return_dict)) return(cluster_colors)
        messager("Adding metadata new column to obj:",color_col)
        obj@meta.data[[color_col]] <- cluster_colors[obj@meta.data[[cluster_col]]]
        return(obj)
    }
}