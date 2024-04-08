#' Get variance explained
#' 
#' Get the proportion of total variance explained
#'  by a given dimensionality reduction 
#'  (0-1 where 1 indicates 100% of variance explained). 
#' @export
#' @examples
#' obj <- get_HPO()
#' out <- get_variance_explained(obj)
get_variance_explained <- function(obj,
                                   reduction = names(obj@reductions)[1],
                                   layer = "scale.data",
                                   dims=NULL
                                   ){
    dr <- obj[[reduction]] 
    mat <- Seurat::GetAssayData(obj, layer = layer)
    total_variance <- sum(matrixStats::rowVars(mat)) 
    ## EigenValues
    if(length(dr@stdev)==0){
        messager("Computing stdev for reduction.")
        dr@stdev <- apply(obj[[reduction]]@cell.embeddings,2,
                          stats::sd)|>unname()
    }
    eigValues <- (dr@stdev)^2|>`names<-`(colnames(dr@cell.embeddings)) 
    if(!is.null(dims)) {
        if(length(dims)>length(eigValues) || 
           max(dims)>length(eigValues) ){
            stopper(
                "The length of dims must be less than or equal to",
                "the number of dimensions."
                )
        }  
        eigValues <- eigValues[dims]
    }
    varExplained <- eigValues / total_variance
    messager("Proportion of total variance explained by",
             length(eigValues),"dimensions:",
             round(sum(varExplained),3))
    return(varExplained)
}