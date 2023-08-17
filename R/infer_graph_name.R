infer_graph_name <- function(obj, 
                             graph_name = NULL,
                             assay = NULL,
                             keys = NULL,
                             ignore_has_graph = TRUE,
                             verbose = TRUE){
    if(is.null(graph_name)){
        if(!is.null(keys)) {
            graph_name <- paste(keys,"cor",sep="_") 
        } else if (!is.null(assay)){
            graph_name <- paste(assay,"cor",sep="_") 
        } else {
            assay <- Seurat::Assays(obj)[1]
            graph_name <- paste(assay,"cor",sep="_") 
        }
    }
    graph_name <- graph_name[1]
    #### Don't check whether graph exists. Always return some graph_name ####
    if(ignore_has_graph){
        messager("Inferred graph_name:",graph_name,v=verbose)
        return(graph_name)
    }
    #### Check whether graph exists. If not, return NULL ####
    if(has_graph(obj = obj, 
                 graph_names = graph_name)){
        return(graph_name)
    } else if (has_graph(obj = obj, 
                         graph_names = "cor")) {
        return("cor")
    } else {
        return(NULL)
    } 
}
