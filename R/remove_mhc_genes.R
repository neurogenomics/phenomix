remove_mhc_genes <- function(dat,
                             verbose = TRUE,
                             ...) {
    MHC_genes <- get_mhc_genes(
        verbose = verbose,
        ...
    )
    MHC_hits <- data.table::uniqueN(dat$GENE[!dat$GENE %in% unique(MHC_genes$SYMBOL)])
    messager("Removing", length(MHC_hits), "MHC gene(s) from data.",
        v = verbose
    )
    dat <- dat[!GENE %in% MHC_hits, ]
    return(dat)
}
