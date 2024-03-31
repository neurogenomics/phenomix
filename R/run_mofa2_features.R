run_mofa2_features <- function(obj,
                               features,
                               assay){
    if(length(assay)>1){ 
        reduce <- NULL
    } else {
        reduce <- union
    }
    if(is.null(features)) {
        features <- rownames(obj) 
    } else if(all(features=="variable_features")){
        features <- scKirby::get_variable_features(obj = obj, 
                                                   reduce = reduce)
    }
    ## Ensure you duplicate the list to match the number of Seurat assays
    if(length(assay)>1 &&
       !is.list(features)){ 
        features <- lapply(obj@assays,function(...)features)
    }   
    return(features)
}