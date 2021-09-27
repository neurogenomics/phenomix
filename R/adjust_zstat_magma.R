#' Adjust GWAS Z-statistic: MAGMA
#'
#' Adjust MAGMA Z-statistic from \emph{.genes.out} files
#'
#' Used when you want to directly analyse the gene-level
#' Z-scores corrected for gene length,  etc.
#'
#' @param magma MAGMA \emph{.genes.out} data.table
#' produced by \link[MAGMA.Celltyping]{map.snps.to.genes}.
#' @param drop_MHC Drop genes from the MHC
#' (Major Histocompatibility Complex) region.
#' @param model Statistical model to use.
#' Defaults to \link[stats]{lm}.
#' @param formula Formula to use in \code{model}.
#' @param ... Additional arguments passed to \code{model}.
#' @inheritParams stats::p.adjust
#'
#' @keywords internal
#' @importFrom stats lm p.adjust
#' @importFrom data.table setkey :=
adjust_zstat_magma <- function(magma,
                               drop_MHC = TRUE,
                               method = "bonferroni",
                               model = NULL,
                               formula = ZSTAT ~ NSNPS + logNSNPS + NPARAM +
                                   logNPARAM + GENELEN + logGENELEN,
                               verbose = TRUE,
                               ...) {
    if (is.null(model)) {
        messager("Defaulting to model stats::lm.",
            v = verbose
        )
        model <- stats::lm
    }
    data.table::setkey(magma, "P")
    magma[, Q := stats::p.adjust(P, method = method)]
    magma[, logNSNPS := log(NSNPS)]
    magma[, logNPARAM := log(NPARAM)]
    #### Drop MHC genes ####
    if (drop_MHC) {
        message("Dropping genes from the MHC region.")
        # DROP MHC: chr6, 25-34 MB
        # Which genome build?
        magma <- magma[!(magma$CHR == 6 & magma$START >= 25000000 & magma$STOP <= 34000000), ]
    }
    #### Calculate gene length ####
    magma[, GENELEN := abs(STOP - START)]
    magma[, logGENELEN := log(GENELEN)]

    #### Regress out effects of NSNPS and NPARAM ####
    # (see 'boxplots_by_decile.r' and the section on downsampling for info)
    #--- NSNPS only really has mjaor effects (i.e. zscore+2) when a gene has ~10000 SNPS
    message("Adjusting MAGMA ZSTAT.")
    mod <- model(
        formula = formula,
        data = magma,
        ...
    )
    magma$ADJ_ZSTAT <- magma$ZSTAT - (magma$NSNPS * mod$coefficients[2] + magma$logNSNPS * mod$coefficients[3] +
        magma$NPARAM * mod$coefficients[4] + magma$logNPARAM * mod$coefficients[5] +
        magma$GENELEN * mod$coefficients[6] + magma$logGENELEN * mod$coefficients[7]
    )
    #### Not the same as above method! ####
    # magma$predict <- stats::predict(object = mod, data=magma$ZSTAT)
    # qplot( magma$ADJ_ZSTAT, magma$predict, alpha=.1) +theme_bw()
    return(magma)
}
