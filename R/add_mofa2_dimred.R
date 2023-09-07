add_mofa2_dimred <- function(obj,
                             model,
                             assay = NULL,
                             keys = NULL,
                             verbose = TRUE) {
    # devoptera::args2vars(add_mofa2_dimred)
    
    DRL <- scKirby::mofa_to_dro(obj = model,
                                keys = keys,
                                verbose = verbose)
    was_list <- if(!is.list(obj)){
        obj <- list(obj)
        FALSE 
    } else {
        TRUE
    }
    obj <- lapply(obj, function(o){
        if(scKirby::is_class(obj,"seurat")){
            for(d in names(DRL)){
                o[[d]] <- DRL[[d]]
            }
        } 
        o
    })
    #### Return ####
    if(isFALSE(was_list)) obj <- obj[[1]] 
    return(obj)
}
