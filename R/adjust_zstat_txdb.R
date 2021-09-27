#' Adjust GWAS Z-statistic: txdb
#'
#' Adjust the GWAS Z-statistic for confounders like gene length.
#'
#' Used when you want to directly analyse the gene-level
#' Z-scores corrected for gene length,  etc.
#'
#' @param dat \code{data.table}
#' produced by \link[phenomix]{map_snps2genes_txdb}.
#' @param drop_MHC Drop genes from the MHC
#' (Major Histocompatibility Complex) region.
#' @param model Statistical model to use.
#' Defaults to \link[stats]{lm}.
#' @param formula Formula to use in \code{model}.
#' @param ... Additional arguments passed to \code{model}.
#' @inheritParams stats::p.adjust
#'
#' @return \code{dat} with the new column "ADJ_ZSTAT".
#'
#' @export
#' @importFrom stats lm p.adjust
#' @importFrom data.table setkey
adjust_zstat_txdb <- function(dat,
                              drop_MHC = TRUE,
                              method = "bonferroni",
                              model = NULL,
                              formula = ZSTAT ~ NSNPS + logNSNPS +
                                  GENELEN + logGENELEN,
                              verbose = TRUE,
                              ...) {
    if (is.null(model)) {
        messager("Defaulting to model stats::lm.",
            v = verbose
        )
        model <- stats::lm
    }

    data.table::setkey(dat, "P_mean")
    dat[, Q := stats::p.adjust(P_mean, method = method)]
    dat[, logNSNPS := log(NSNPS)]
    # dat[,logNPARAM:=log(NPARAM)]
    #### Drop MHC genes ####
    if (drop_MHC) {
        message("Dropping genes from the MHC region.")
        # DROP MHC: chr6, 25-34 MB
        MHC_genes <- get_mhc_genes()
        dat <- dat[!SYMBOL %in% unique(MHC_genes$SYMBOL), ]
    }
    #### Calculate gene length ####
    dat[, logGENELEN := log(GENELEN)]

    #### Regress out effects of NSNPS and NPARAM ####
    # (see 'boxplots_by_decile.r' and the section on downsampling for info)
    #--- NSNPS only really has mjaor effects (i.e. zscore+2) when a gene has ~10000 SNPS
    messager("Adjusting dat ZSTAT.", v = verbose)
    mod <- model(
        formula = formula,
        data = dat,
        ...
    )
    dat$ADJ_ZSTAT <- dat$ZSTAT - (dat$NSNPS * mod$coefficients[2] +
        dat$logNSNPS * mod$coefficients[3] +
        dat$GENELEN * mod$coefficients[4] +
        dat$logGENELEN * mod$coefficients[5]
    )
    #### Not the same as above method! ####
    # dat$predict <- stats::predict(object = mod, data=dat$ZSTAT)
    # qplot( dat$ADJ_ZSTAT, dat$predict, alpha=.1) +theme_bw()
    return(dat)
}
