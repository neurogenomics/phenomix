phenomix_query_magma <- function(ids = NULL,
                                 meta = opengwas_meta(ids = ids),
                                 as_matrix = TRUE, 
                                 verbose = TRUE,
                                 ...){
    magma_out <- NULL;
    
    # meta <- meta[grepl("eqtl-a-",id),][seq(10)]
    meta <- meta[RCurl::url.exists(magma_out),]
    magma_out <- stats::setNames(meta$magma_out,meta$id)
    res <- magma_matrix(magma_out = magma_out,
                      as_matrix = as_matrix,
                      verbose = verbose,
                      ...)  
    return(res)
}