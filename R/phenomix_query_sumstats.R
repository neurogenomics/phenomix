phenomix_query_sumstats <- function(ids = NULL,
                                    meta = opengwas_meta(ids = ids),
                                    save_dir = tempdir(),
                                    query_granges = NULL, 
                                    as_matrix = FALSE, 
                                    verbose = TRUE,
                                    ...){
    MAP <- function(range, 
                    file, 
                    save_dir=tempdir(),
                    verbose=TRUE, 
                    ...){
        id <- gsub("\\.tsv|\\.bgz","",basename(file))
        messager("Querying:",id,parallel=TRUE)
        if(!RCurl::url.exists(file)){
            messager("WARNING: URL for",id,"does not exist. Skipping.",
                     v=verbose)
            return(NULL)
        }                              
        if(is.null(range)){
            echotabix::read_bgz(path = file,
                                verbose = verbose,
                                ...)
        } else {
            echotabix::query(target_path = file, 
                             query_granges = range, 
                             query_save_path = file.path(
                                 save_dir,paste0(id,".tsv.gz")
                             ),
                             verbose = verbose,
                             ...)
        } 
    }
    #### Query ####
    files <- stats::setNames(meta$url,meta$id)
    res <- GenomicFiles::reduceFiles(
        ranges = query_granges,
        files = files,
        MAP = MAP,
        REDUCE = function(mapped,...){data.table::rbindlist(mapped)},
        save_dir = save_dir,
        verbose = verbose)
    # res <- lapply(X = files,
    #               FUN = MAP,
    #               range = query_granges,
    #               save_dir = save_dir,
    #               verbose = verbose)
    closeAllConnections() 
    #### Merge and return ####
    if(isTRUE(as_matrix)){
        return(
            phenomix_merge(res = res,
                           as_matrix = as_matrix,
                           verbose = verbose)
        )
    } else {
        return(res) 
    } 
}