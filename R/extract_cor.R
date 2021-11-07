#' Extract correlation matrix
#' 
#' Extract correlation matrix stored as a graph.
#' 
#' @param assay Assay to used to compute correlation graph 
#' if one does not already exist.
#' @param graph_name Name of the graph to extract.
#' @param method Pairwise correlation method.
#' @inheritParams extract_embeddings
#' @inheritParams SeuratObject::Reductions
#' 
#' @return Trait-trait correlation matrix.
#' 
#' @export
extract_cor <- function(obj,
                        assay = NULL,
                        slot = NULL,
                        reduction = NULL,
                        graph_name = NULL,
                        method = "pearson",
                        verbose = TRUE) {
    if (is_seurat(obj = obj)) {
        #### Search for (cor) graphs ####
        graphs <- names(obj@graphs)
        if (is.null(graph_name)) {
            if (length(graphs) > 0) {
                cor_graphs <- grep("*_cor", graphs)
                if (length(cor_graphs) > 0) {
                    graph_name <- cor_graphs[1]
                } else {
                    cmat <- compute_cor(
                        obj = obj,
                        assay = assay,
                        slot = slot,
                        reduction = reduction,
                        method = method,
                        return_obj = FALSE
                    )
                    return(cmat)
                }
            } else {
                #### If no graphs, compute cor graph ####
                cmat <- compute_cor(
                    obj = obj,
                    assay = assay,
                    slot = slot,
                    reduction = reduction,
                    return_obj = FALSE
                )
                return(cmat)
            }
        }
        messager("Using graph:", graph_name, v = verbose)
        cmat <- extract_graph(
            obj = obj,
            graph_name = graph_name,
            verbose = verbose
        )
        return(cmat)
    } else {
        cmat <- compute_cor(
            obj = obj,
            assay = assay,
            slot = slot,
            reduction = reduction,
            return_obj = FALSE
        )
        return(cmat)
    }
}
