#' Merge phenomix queries
#' 
#' Merge results from \link[phenomix]{query_phenomix} 
#' into a single sparse matrix or \link[data.table]{data.table} filled with summary statistics.
#' @param res Results from \link[phenomix]{query_phenomix}. 
#' Can be a list of \link[data.table]{data.table}s or sparse matrices.
#' @param value_var Column to fill the matrix with.
#' @param key_var A column (or set of columns) to be used when
#'  merging datasets together.
#' @param save_path Path to save matrix to in RDS format
#' (set to \code{NULL} to skip this step).
#' @param as_matrix Convert the data to a sparse matrix.
#' @param verbose Print messages.
#' @inheritParams scKirby::to_sparse
#' @inheritParams data.table::merge.data.table
#' @returns Sparse matrix
#' 
#' @export
#' @import scKirby
#' @import data.table
#' @examples
#' query_granges <- GenomicRanges::GRanges(c("2:15000-17000","3:60000-70000"))
#' ids <- tail(phenomix::OpenGWAS[!grepl("^finn",id),]$id,3)
#' res <- phenomix_query(ids=ids, query_granges=query_granges)
#' X <- phenomix_merge(res) 
phenomix_merge <- function(res,
                           value_var="BETA",
                           key_var=c("SNP","A1","A2","CHR","BP"),
                           save_path=NULL,
                           as_matrix=TRUE,
                           all=TRUE,
                           verbose=TRUE){
    
    #### Remove NULLs ####
    res <- rm_empty(res)
    #### Sparse matrix: merge ####
    if(scKirby::is_class(res[[1]],"sparse_matrix")){
        X <- Reduce(Seurat::RowMergeSparseMatrices,res)
    } else {
    #### data.table: Subset and merge ####
        X <- lapply(names(res), function(nm){ 
            res[[nm]][,key_:= do.call(paste, c(.SD, sep = "_")),
                      .SDcols = key_var][,c("key_",value_var),with=FALSE,
                                         keyby="key_"] |>
                data.table::setnames(value_var,nm)
        }) |> 
            Reduce(f = function(x, y) merge(x, y, all = all)) 
        #### Convert to sparse matrix ####
        if(isTRUE(as_matrix)){
            X <- scKirby::to_sparse(obj = X, 
                                    verbose = verbose) 
        }
    } 
    #### Report ####
    messager("Data dimensions:",paste(dim(X),collapse = " x "),v=verbose)
    #### Save ####
    if(!is.null(save_path)){
        messager("Saving merged data -->",save_path,v=verbose)
        dir.create(dirname(save_path),showWarnings = FALSE, recursive = TRUE)
        saveRDS(X,save_path)
    }
    return(X)
}
