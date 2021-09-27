#' Aggregate summary statistics
#'
#' Efficiently aggregate GWAS/QTL summary statistics
#' according to some grouping column (e.g."GENEID").
#'
#' Requires the columns
#'
#' @param merged_hits \link[data.table]{data.table}
#' with the columns "BETA", "P","SE", and \code{agg_var}.
#' @param drop_na Drop rows with \code{NA} in any column.
#' @param agg_var Variable to aggregate \code{merged_hits} by.
#' @param verbose Print messages.
#'
#' @return \link[data.table]{data.table} with
#' mean and standard deviation of "BETA", "P", and "SE", as well as
#' number of rows and number of unique SNPs per grouping level.
#'
#' @source \href{https://stackoverflow.com/questions/59059743/apply-different-functions-to-different-columns-programmatically-in-data-table-r}{data.table aggregation solution}
#'
#' @export
#' @importFrom data.table uniqueN
#' @importFrom stats complete.cases
aggregate_sumstats <- function(merged_hits,
                               drop_na = TRUE,
                               agg_var = "GENEID",
                               extra_cols = c("ABC.Score", "CellType"),
                               sort_rows = TRUE,
                               verbose = TRUE) {

    #### Check grouping column exists ####
    if (!agg_var %in% colnames(merged_hits)) {
        stop("agg_var must be a column name in merged_hits")
    }
    messager("Aggregating merged_hits by", paste0(agg_var, "."), v = verbose)
    #### Define functions ####
    prepare_agg_funcs_out <- prepare_agg_funcs(
        merged_hits = merged_hits,
        extra_cols = extra_cols
    )
    agg_functions <- prepare_agg_funcs_out$agg_functions
    col_selection <- prepare_agg_funcs_out$col_selection
    #### Aggregate efficiently ####
    gene_hits <- merged_hits[, mapply(
        function(f, x) as.list(f(x)),
        agg_functions, .SD
    ),
    agg_var,
    .SDcols = col_selection
    ]
    #### Sort ####
    if (sort_rows) {
        messager("Sorting rows alphanumerically by", paste0(agg_var, "."),
            v = verbose
        )
        data.table::setkeyv(gene_hits, agg_var)
    }
    ##### Drop NAs ####
    incomplete <- sum(!stats::complete.cases(gene_hits))
    if (drop_na & incomplete > 0) {
        messager("Dropping", formatC(incomplete, big.mark = ","),
            "rows with NAs in any column.",
            v = verbose
        )
        gene_hits <- na.omit(gene_hits)
    }
    #### Report ####
    messager("Returning gene_hits data.table with",
        formatC(nrow(gene_hits), big.mark = ","), "rows.",
        v = verbose
    )
    return(gene_hits)
}



#
#
# gene_hits <- merged_hits[,list(BETA_mean=mean(BETA,na.rm=TRUE),
#                                # SE_mean=mean(SE,na.rm=TRUE),
#                                P_mean=mean(P,na.rm=TRUE),
#
#                                BETA_sd=sd(BETA,na.rm=TRUE),
#                                SE_sd=sd(SE,na.rm=TRUE),
#                                P_sd=sd(P,na.rm=TRUE),
#                                P_sd=sd(P,na.rm=TRUE),
#
#                                # regions=data.table::uniqueN(region, na.rm=TRUE),
#                                nonuniqueSNPS=.N,
#                                NSNPS= data.table::uniqueN(SNP, na.rm=TRUE)
# ),
# by=c(agg_var)]
