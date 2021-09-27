extract_graph <- function(obj,
                          graph_name=NULL,
                          verbose=TRUE){
    if(methods::is(obj,"Seurat")){
        if(is.null(graph_name)) graph_name <- rev(names(obj@graphs))[1] 
        if(any(graph_name %in% names(obj@graphs))){
            messager("Using graph:",graph_name,v=verbose)
        }
        fgraph <- obj@graphs[[graph_name]]
        return(fgraph)
    } else if (methods::is(obj,"Graph")) {
        messager("Using obj as graph.",v=verbose)
        fgraph <- obj
        
    } else {
        messager("No graph found. Returning NULL.",v=verbose)
        fgraph <- NULL
    }
    return(fgraph)
}
