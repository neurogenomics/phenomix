#' Construct gene x trait matrix from MAGMA results
#'
#' Once you've run \link[MAGMA.Celltyping]{map_snps_to_genes} 
#' on multiple GWAS/QTL,
#' you can aggregate the gene-wise results into a single gene x trait matrix.
#'
#' @param magma_out Paths to \emph{genes.out} files
#'  produced by \link[MAGMA.Celltyping]{map_snps_to_genes}.
#' @param metric Which metric to fill the matrix with
#' \itemize{
#'  \item{\code{"ADJ_ZSTAT"} : }{MAGMA Z-statistic after adjusting for
#'  gene length and other confounds (\emph{DEFAULT}).}
#'  \item{\code{"ZSTAT"} : }{Unadjusted MAGMA Z-statistic.}
#'  \item{\code{"P"} : }{MAGMA P-value.}
#' }
#' @param nCores Number of cores to parallelise across.
#' @param fillna Value to fill \code{NA}s with.
#' @param trans_fun Transformation function. 
#' @param agg_fun Aggregation function.
#' @param save_path Path to save results as RDS object.
#' @param verbose Print messages.
#' @inheritParams adjust_zstat
#' @inheritParams phenomix_merge
#' @return Sparse matrix.
#' 
#' @export 
#' @import data.table
#' @examples
#' meta <- opengwas_meta()
#' meta <- meta[grepl("eqtl-a-",id),][seq(5)]
#' X <- magma_matrix(magma_out = meta$magma_out, workers=10) 
magma_matrix <- function(magma_out,
                         metric = c("ADJ_ZSTAT", "ZSTAT", "P"),
                         drop_mhc = TRUE,
                         save_path = NULL, 
                         fillna = FALSE,
                         trans_fun = NULL,#function(x){-log10(x)},
                         agg_fun = "mean",
                         all = TRUE,
                         as_matrix = TRUE,
                         workers = 1,
                         verbose = TRUE) {
    # devoptera::args2vars(magma_matrix) 
    
    GENE <- SYMBOL <- NULL;
    metric <- toupper(metric[1])
    if(is.null(names(magma_out))){
        names(magma_out) <- basename(dirname(dirname(dirname(magma_out))))
    }
    fill_value <- if (metric == "P") 1 else 0
    #### Merge all results ####
    BPPARAM <- assign_cores(workers = workers)
    DAT <- BiocParallel::bplapply(BPPARAM = BPPARAM,
        X = seq(length(magma_out)),
        FUN = function(i) {
            x <- names(magma_out)[i]
            tryCatch({
                messager(x," : ",i," / ",length(magma_out), parallel=TRUE,
                         v= workers==1)
                dat <- data.table::fread(magma_out[[x]])
                #### Adjust ZSTAT ####
                if (metric == "ADJ_ZSTAT") {
                    dat <- adjust_zstat(
                        dat = dat,
                        drop_mhc = drop_mhc,
                        gene_col = "GENE",
                        log_vars = c("NSNPS", "NPARAM", "GENELEN"),
                        formula = ZSTAT ~ NSNPS + logNSNPS + NPARAM +
                            logNPARAM + GENELEN + logGENELEN,
                        verbose = workers==1
                    )
                }
                dat[,GENE:=as.character(GENE)]
                data.table::setkey(dat, "GENE")
                data.table::setnames(dat, old = metric, new = x)
                dat[, c("GENE",x), with=FALSE]
            }, error=function(e){messager(e,v=verbose);NULL})
        }
    )
    DAT <- rm_empty(res = DAT)
    DAT <- DAT |> Reduce(f = function(x, y)  merge(x, y, all=all))
    #### Makes colnames compatible with data.table and Matrix format ####
    messager(
        "Replacing '-'/':'/'.' with '_' in colnames ",
        "to make compatible with sparse matrix format",
        v=verbose
    )
    old_colnames <- colnames(DAT)
    new_colnames <- gsub("-|[:]|[.]", "_", colnames(DAT))
    data.table::setnames(DAT, old_colnames, new_colnames)
    #### Fill NAs ####
    if (isTRUE(fillna)) {
        messager("Replacing NAs with",fill_value,".",v=verbose)
        #### Very efficient NA replacement  ####
        na.replace <- function(v, value = fill_value) {
            v[is.na(v)] <- value
            v
        }
        for (i in names(DAT)) {
            eval(parse(text = paste("DT[,", i, ":=na.replace(", i, ")]")))
        }
    }
    #### Translate gene IDs ####    
    DAT <- translate_geneids_txdb(
        gene_hits = DAT, 
        gene_col = "GENE",
        verbose = verbose
    )
    DAT <- DAT[,GENE:=SYMBOL][,-c("SYMBOL")]
    #### Return early as data.table ####
    if(isFALSE(as_matrix)) return(DAT)
    #### Convert to sparse matrix ####
    X <- scKirby::to_sparse(DAT, verbose = verbose)
    #### -log transform ####
    if (!is.null(trans_fun)) {
        if(isTRUE(verbose)){ 
            messager("Applying transformation to matrix:", 
                     v=verbose) 
            methods::show(trans_fun)
        }
        X <- trans_fun(X)
    } 
    if (!is.null(agg_fun) &&
        (sum(duplicated(rownames(X))) > 0)) {
        messager("Aggregating duplicated gene rows by ", agg_fun,v=verbose) 
        X <- orthogene:::aggregate_rows(
            X = X,
            groupings = rownames(X),
            agg_fun = agg_fun,
            as_DelayedArray = FALSE
        )
        messager("Aggregated matrix dimensions:",
                 paste(formatC(dim(magma_matrix),big.mark = ","),
                       collapse = " x "),v=verbose)
    }
    #### Save ####
    if (!is.null(save_path)) {
        messager("Saving results -->", save_path,v=verbose)
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        saveRDS(X, save_path)
    }
    return(X)
}
