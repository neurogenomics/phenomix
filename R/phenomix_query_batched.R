#' Query \href{https://phenomix.dsi.ic.ac.uk/}{phenomix}
#' 
#' Query the \href{https://phenomix.dsi.ic.ac.uk/}{phenomix} database
#'  for harmonised GWAS or QTL datasets. 
#' @param batch_size Number of datasets to process per batch.
#' @inheritParams opengwas_meta
#' @inheritParams phenomix_merge
#' @inheritParams phenomix_query
#' @inheritParams MungeSumstats::find_sumstats
#' @inheritParams echotabix::query_table
#' @inheritParams data.table::merge.data.table
#' @inheritDotParams echotabix::query_table 
#' @returns data.table of genomic regions
#' 
#' @export 
#' @importFrom stats setNames
#' @examples
#' query_granges <- GenomicRanges::GRanges(c("2:15000-17000","3:60000-70000"))
#' ids <- tail(phenomix::OpenGWAS[!grepl("^finn",id),]$id,6)
#' X <- phenomix_query_batched(ids=ids,
#'                             query_granges=query_granges,
#'                             batch_size=2)
phenomix_query_batched <- function(ids = NULL,
                                   meta = opengwas_meta(ids = ids), 
                                   save_dir = tempdir(),
                                   query_granges = NULL, 
                                   batch_size=100,
                                   as_matrix=TRUE,
                                   value_var="BETA",
                                   key_var=c("SNP","A1","A2","CHR","BP"),
                                   all=TRUE, 
                                   verbose=TRUE,
                                   ...){
    
    # devoptera::args2vars(phenomix_query_batched)
     
    batches <- split(meta$id, ceiling(seq_along(meta$id)/batch_size))
    tmp <- tempfile()
    dtl <- lapply(seq(length(batches)), function(i){
        messager("Processing batch: ",i,"/",length(batches),
                 parallel = TRUE, v = verbose)
        save_path <- paste0(tmp,"_phenomatrix_batch",i,".rds")
        phenomix_query(meta = opengwas_meta(ids = batches[[i]]), 
                       save_dir = save_dir,
                       query_granges = query_granges, 
                       query_save = FALSE, 
                       verbose = verbose) |>
            phenomix_merge(value_var = value_var,
                           key_var = key_var,
                           save_path = NULL,
                           as_matrix = FALSE,
                           verbose = verbose)
    })
    #### Merge data.tables #### 
    messager("Merging",length(dtl),"data.tables.",v=verbose)
    X <- Reduce(f = function(x, y) merge(x, y, all = all), x = dtl) 
    #### Convert to sparse matrix ####
    if(isTRUE(as_matrix)){
        X <- scKirby::to_sparse(obj = X,
                                verbose = verbose)
    }
    return(X)
}
