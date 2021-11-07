#' Extract correlation matrix
#' 
#' Extract correlation matrix stored as a graph.
#' 
#' @param assay Assay to used to compute correlation graph 
#' if one does not already exist.
#' @param graph_name Name of the graph to extract.
#' @param method Pairwise correlation method.
#' @inheritParams compute_cor
#' @inheritParams extract_embeddings
#' @inheritParams SeuratObject::Reductions
#' 
#' @return Trait-trait correlation matrix.
#' 
#' @export
#' hpo <- phenomix::get_HPO()
#' cmat <- extract_cor(obj = hpo,
#'                     reduction = "pca")
extract_cor <- function(obj,
                        assay = NULL,
                        slot = NULL,
                        reduction = NULL,
                        graph_name = NULL,
                        method = "pearson",
                        return_obj = TRUE,
                        verbose = TRUE) {
    
    #### When obj is a Seurat object ####
    if (is_seurat(obj = obj)) {
        #### Reassign name to distinguish from other cor matrices ####
        graph_name <- infer_graph_name(obj = obj, 
                                       graph_name = graph_name,
                                       assay = assay, 
                                       reduction = reduction, 
                                       ignore_has_graph = TRUE)
        #### Check if graph exists ####
        if (has_graph(obj = obj, 
                      graph_names = graph_name)){
            #### Use existing _cor graph #### 
            messager("Using pre-computed graph:", graph_name, v = verbose)
            cmat <- extract_graph(
                obj = obj,
                graph_name = graph_name,
                verbose = verbose) 
        } else {
            #### Compute new _cor graph #### 
            out <- compute_cor(
                obj = obj,
                assay = assay,
                slot = slot,
                reduction = reduction,
                method = method,
                return_obj = return_obj
            ) 
            return(out)
        }  
        #### Return results ####
        if(return_obj){
            messager("Returning Seurat object.",v=verbose) 
            obj@graphs[[graph_name]] <- cmat 
            return(obj)
            
        } else {
            messager("Returning correlation matrix.",v=verbose)
            return(cmat)   
        }
    } else {
        #### When obj is a matrix ####
        if(isTRUE(return_obj)){
            messager("Cannot return object when obj is a matrix.",
                     "Returning correlation matrix instead",v=verbose)
        }
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
