#' Query \href{https://phenomix.dsi.ic.ac.uk/}{phenomix}
#' 
#' Query the \href{https://phenomix.dsi.ic.ac.uk/}{phenomix} database
#'  for harmonised GWAS or QTL datasets.
#' 
#' @param save_dir Directory to save query results to.
#' @param data_type Data type to query:
#' \itemize{
#' \item{"sumstats" : }{Genome-wide SNP-level summary statistics 
#' that have previously been harmonised using 
#' \link[MungeSumstats]{format_sumstats}.}
#' \item{"magma" : }{Genome-wide gene-level 
#' \href{https://ctg.cncr.nl/software/magm}{MAGMA} results.}
#' }
#' @inheritParams opengwas_meta
#' @inheritParams phenomix_merge
#' @inheritParams MungeSumstats::find_sumstats
#' @inheritParams echotabix::query_table
#' @inheritDotParams echotabix::query_table
#' @returns data.table of genomic regions
#' 
#' @export 
#' @import GenomicFiles
#' @importFrom stats setNames
#' @examples
#' query_granges <- GenomicRanges::GRanges(c("2:15000-16000","3:60000-61000"))
#' res <- phenomix_query(ids="bbj-a-1", query_granges=query_granges)
phenomix_query <- function(ids = NULL,
                           meta = opengwas_meta(ids = ids),
                           save_dir = tempdir(),
                           query_granges = NULL,
                           data_type = c("sumstats","magma"),
                           as_matrix = FALSE, 
                           verbose = TRUE,
                           ...){
    # devoptera::args2vars(phenomix_query)
    requireNamespace("echotabix")
     
    data_type <- tolower(data_type)[1]
    if(is.character(query_granges) && 
       query_granges=="cs2g_ukb"){
        query_granges <- get_cs2g(as_granges = TRUE,
                                  verbose = verbose)$data
    }
    #### Query SNP-level GWAS/QTL sumstats ####
    if(data_type=="sumstats"){
        res <- phenomix_query_sumstats(ids = ids,
                                       meta = meta,
                                       save_dir = save_dir,
                                       query_granges = query_granges, 
                                       as_matrix = as_matrix, 
                                       verbose = verbose,
                                       ...)
    } else if (data_type=="magma"){
        res <- phenomix_query_magma(ids = ids,
                                    meta = meta,
                                    as_matrix = as_matrix, 
                                    verbose = verbose,
                                    ...) 
    } 
    return(res)
}
