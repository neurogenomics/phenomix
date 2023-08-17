p_adjust <- function(dat,
                     method = "bonferroni",
                     P_var = NULL,
                     Q_var = "Q",
                     verbose = TRUE) {
    Q <- NULL;
    
    if("Q" %in% names(dat)){
        messager("Using existing Q column.",v=verbose) 
    } else {
        #### Figure out which column to use for P ####
        if (all(!is.null(P_var), P_var %in% colnames(dat))) {
            P_use <- P_var
        } else if ("P" %in% colnames(dat)) {
            P_use <- "P"
        } else if ("p" %in% colnames(dat)) {
            P_use <- "P"
        } else if ("P_mean" %in% colnames(dat)) {
            P_use <- "P_mean"
        } else {
            stopper("No usable P_var detected.") 
        }
        #### Report ####
        messager("Computing", Q_var, "by correcting", P_use, "with",
                 paste0(method, "."),
                 v = verbose
        )
        #### Compute Q ####
        data.table::setkeyv(dat, P_use)
        dat[, Q := stats::p.adjust(get(P_use), method = method)]
        data.table::setnames(dat, Q_var, "Q")
    } 
}
