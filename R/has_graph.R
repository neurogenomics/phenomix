has_graph <- function(obj,
                      graph_names = NULL){
    if(scKirby::is_class(obj,"seurat")){
        #### Check if "graphs" slot even exists ####
        if(!"graphs" %in% methods::slotNames(obj)){
            return(FALSE)
        }
        #### Check if specific graphs exist ####
        all_graphs <- names(obj@graphs)
        if(all(is.null(graph_names))){
            # If no graph_name given, 
            # just tell us whether there's any graphs.
            return(length(all_graphs)>0)
        } else {
            # If graph_name given, tell us whether it exists in the obj
            return(any(graph_names %in% all_graphs))
        }
    } else{ return(FALSE) }
}