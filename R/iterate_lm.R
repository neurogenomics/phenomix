#' Iterate linear regression across all traits and celltypes
#'
#' @param xmat gene x trait matrix.
#' @param ymat gene x celltype matrix.
#' @param correction_method Multiple-testing correction 
#' method to be passed to \code{stats::p.adjust}.
#' @param qvalue_thresh q-value threshold to use when report
#'  significant results summary.
#' @param quantize A named list where the values of "x" and "y" indicaite the 
#' number of quantiles to bin the respective "xmat" and "ymat" datasets into.
#' @inheritParams set_cores
#'
#' @export
#' @import data.table 
#' @importFrom stats p.adjust
#' @examples
#' ### Phenotype dataset
#' obj <- get_HPO()
#' xmat <- scKirby::get_x(obj)[["RNA.data"]]
#' xmat <- xmat[,seq(10)] # Let's use just 10 traits as an example
#'
#' ### Celltypes Dataset
#' ctd <- get_ctd()
#' ymat <- ctd[[1]]$specificity
#' res_lm <- iterate_lm(xmat = xmat,
#'                      ymat = ymat, 
#'                      workers = 1)
iterate_lm <- function(xmat,
                       ymat,
                       correction_method = "BH",
                       qvalue_thresh = .05,
                       quantize = list(x=NULL,
                                       y=NULL),
                       progressbar = TRUE,
                       workers = NULL,
                       verbose=TRUE) {
    # devoptera::args2vars(iterate_lm)
    
    x <- y <- p <- q <- NULL;
    data.table::setDTthreads(threads = 1)
    cores <- set_cores(workers = workers,
                       progressbar = progressbar,
                       verbose = verbose) 
    ## Filter data
    X_list <- filter_matrices(X_list = list(xmat=xmat,
                                            ymat=ymat),
                              verbose = verbose)
    xmat <- X_list$xmat
    ymat <- X_list$ymat 
    ### Run lm for all combinations of xmat and ymat
    messager(
        "Running", formatC(ncol(xmat) * ncol(ymat), big.mark = ","),
        "tests:", formatC(ncol(xmat), big.mark = ","), "xmat columns x",
        formatC(ncol(ymat), big.mark = ","), "ymat columns.",v=verbose
    ) 
    ## Quantize data
    if(is.numeric(quantize$x)){
        xmat <- scKirby:::quantize_matrix(X=xmat,
                                         n=quantize$x, 
                                         verbose = verbose)
    }
    if(is.numeric(quantize$y)){
        ymat <- scKirby:::quantize_matrix(X=ymat,
                                         n=quantize$y, 
                                         verbose = verbose)
    }
    
    ## Iterate over xmat cols
    lm_res <- BiocParallel::bplapply(
        BPPARAM = cores$params,
        X = stats::setNames(seq(ncol(xmat)), 
                            colnames(xmat)), 
        FUN = function(i) {
            tt <- colnames(xmat)[i]
            if(!progressbar){
                messager("-",tt,": (",i,"/",ncol(xmat),")", 
                         parallel=TRUE)    
            } 
            ## Prepare data for rstatix 
            dt <- melt_merge_matrices(
                xmat = xmat[,tt, drop=FALSE],
                ymat = ymat) 
            dt <- dt[!is.na(x) & !is.na(y),]
            if(isFALSE(dt_var_check(dt,"feature",verbose=!progressbar)) ||
               isFALSE(dt_var_check(dt, "x", verbose=!progressbar)) ||
               isFALSE(dt_var_check(dt, "y", verbose=!progressbar)) ){
                return(NULL)
            }
            ## Run tests
            res <- dt |>
                rstatix::group_by(yvar) |> 
                rstatix::anova_test(formula = y ~ x) |>
                data.table::data.table() 
            return(res)
    }) |> 
        data.table::rbindlist(idcol = "xvar",
                              fill = TRUE) 
    ### Multiple-testing correction 
    lm_res[,q:=stats::p.adjust(p = p, method = correction_method)] 
    ### Filter only sig results 
    messager(formatC(nrow(lm_res[q<qvalue_thresh,]), big.mark = ","),
            "significant results @",
            correction_method, "<", qvalue_thresh)
    ### Return FULL results (not just sig)
    return(lm_res)
}
