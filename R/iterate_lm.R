#' Iterate linear regression across all traits and celltypes
#'
#' @param xmat gene x trait matrix.
#' @param ymat gene x celltype matrix.
#' @param test_method Association testing method to use.
#' @param correction_method Multiple-testing correction 
#' method to be passed to \link[stats]{p.adjust}.
#' @param qvalue_thresh q-value threshold to use when report
#'  significant results summary.
#' @param quantize A named list where the values of "x" and "y" indicate the 
#' number of quantiles to bin the respective \code{xmat} and 
#' \code{ymat} datasets into.
#' @param multivariate If \code{TRUE}, runs tests with each column in 
#' \code{xmat} as a multivariate predictor in a single model 
#' (one model per column in \code{yvar}). If\code{FALSE},
#' runs tests with each column in \code{xmat} as a univariate predictor in
#' separate models (one model per column in \code{yvar}). 
#' @param ... Additional parameters passed to the statistical test function.
#' @inheritParams set_cores
#'
#' @export
#' @import data.table 
#' @importFrom stats p.adjust
#' @examples
#' ### Celltypes Dataset
#' ctd <- get_ctd()
#' xmat <- ctd[[1]]$specificity
#' 
#' ### Phenotype dataset
#' obj <- get_HPO()
#' ymat <- scKirby::get_x(obj)[["RNA.data"]]
#' ymat <- ymat[,seq(10)] # Let's use just 10 traits as an example
#' 
#' res_lm <- iterate_lm(xmat = xmat,
#'                      ymat = ymat, 
#'                      workers = 1)
iterate_lm <- function(xmat,
                       ymat,
                       test_method = c("glm","anova","lm.ridge","rlm"), 
                       multivariate = FALSE,
                       correction_method = "BH",
                       qvalue_thresh = .05,
                       quantize = list(x=NULL,
                                       y=NULL),
                       progressbar = TRUE,
                       workers = NULL,
                       verbose = TRUE,
                       ...) {  
    p <- q <- NULL;
    test_method <- match.arg(test_method)
    data.table::setDTthreads(threads = 1)
    cores <- set_cores(workers = workers,
                       progressbar = progressbar,
                       verbose = verbose) 
    t1 <- Sys.time()
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
    #### Run tests ####
    lm_res <- iterate_lm_long(xmat = xmat,
                              ymat = ymat, 
                              cores = cores,
                              method = test_method,
                              multivariate = multivariate,
                              ...)
    ### Multiple-testing correction 
    if(test_method=="anova"){  
        lm_res[,q:=stats::p.adjust(p = p, method = correction_method)] 
    } else { 
        grep_terms <- if(isTRUE(multivariate)) "^x[:]" else "^x$"
        model_p <- model_q <- term <- xvar <- NULL;
        lm_res|>data.table::setnames("p.value","p")
        model_res <- lm_res[term %in% c("(Intercept)")]|>
            data.table:::setnames(c("p","estimate","statistic"),
                                  c("model_p","model_estimate","model_statistic"))
        model_res[,model_q:=stats::p.adjust(p = model_p, 
                                            method = correction_method)]  
        lm_res <- merge(lm_res[grepl(grep_terms,term)] ,
                        model_res[,c("model_id","model_p","model_q",
                                     "model_estimate","model_statistic")], 
                        by="model_id")
        if(isTRUE(multivariate)){
            lm_res[,xvar:=gsub("^x\\:xvar","",term)]|>
                data.table::setcolorder("xvar",3)
            ## by=.I is EXTREMELY important! otherwise, the p-values will be 
            ## the minimum across ALL rows.
            lm_res[,q:=ifelse(model_q<0.05,p,min(1,model_q+p)), by=.I] 
        } else {
            lm_res[,q:=stats::p.adjust(p = p, method = correction_method)] 
        }
    }   
    #### Report ####
    messager(formatC(nrow(lm_res[q<qvalue_thresh,]), big.mark = ","),
             "significant results @",
             correction_method, "<", qvalue_thresh)
    methods::show(Sys.time()-t1)
    ### Return FULL results (not just sig)
    return(lm_res)
}
