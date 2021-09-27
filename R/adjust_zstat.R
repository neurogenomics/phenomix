#' Adjust GWAS Z-statistic
#'
#' Adjust the GWAS Z-statistic.
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
#' @param log_vars Variables to perform natural log transformation on first.
#' Only run on variables available in \code{dat}.
#' @param formula Formula to use in \code{model}.
#' @param ... Additional arguments passed to \code{model}.
#' @inheritParams stats::p.adjust
#'
#' @return \code{dat} with the new column "ADJ_ZSTAT".
#'
#' @export
#' @importFrom stats lm p.adjust na.omit
#' @importFrom data.table setkey
adjust_zstat <- function(dat,
                         drop_MHC = TRUE,
                         method = "bonferroni",
                         model = NULL,
                         log_vars = c("NSNPS", "NPARAM", "GENELEN"),
                         formula = ZSTAT ~ NSNPS + logNSNPS + NPARAM +
                             logNPARAM + GENELEN + logGENELEN,
                         verbose = TRUE,
                         ...) {

    # dat <- gene_hits;drop_MHC=TRUE; method="bonferroni";   model=NULL; log_vars=c("NSNPS","NPARAM","GENELEN");    formula=ZSTAT ~ NSNPS + logNSNPS + GENELEN + logGENELEN

    #### Adjust p-value ####
    p_adjust(
        dat = dat,
        method = method,
        P_var = NULL,
        Q_var = "Q"
    )
    #### Log relevant variables ###
    log_transform(
        dat = dat,
        log_vars = log_vars,
        verbose = verbose
    )
    #### Drop MHC genes ####
    if (drop_MHC) {
        dat <- remove_mhc_genes(
            dat = dat,
            verbose = verbose
        )
    }

    dat_og <- dat
    dat <- stats::na.omit(dat, "ZSTAT")
    if (nrow(dat) < 10) {
        messager("WARNING: <10 rows remainings in gene-level data. Unable to compute ADJ_ZSTAT.")
        return(dat_og)
    }

    #### Adjust formula ####
    # Make sure all variables are available in dat #
    formula1 <- fix_formula1(
        formula = formula,
        dat = dat
    )

    #### Regress out effects of NSNPS and NPARAM ####
    # (see 'boxplots_by_decile.r' and the section on downsampling for info)
    #--- NSNPS only really has mjaor effects (i.e. zscore+2) when a gene has ~10000 SNPS

    #### Compute ADJ_ZSTAT (done inplace with data.table) ####
    adjust_zstat_run(
        dat = dat,
        model = model,
        formula = formula1,
        verbose = verbose,
        ...
    )

    # #### Method 21 (original from MAGMA.Celltyping) ####
    # dat$ADJ_ZSTAT2 = dat$ZSTAT - (dat$NSNPS*mod$coefficients[2] +
    #                                  dat$logNSNPS*mod$coefficients[3] +
    #                                  dat$GENELEN*mod$coefficients[4] +
    #                                  dat$logGENELEN*mod$coefficients[5]
    #                              )
    #### Method3 (NOTE:  Not the same as the other two methods!) ####
    # dat$ADJ_ZSTAT3 <- stats::predict(object = mod, data=dat$ZSTAT)
    # #### Compare ####
    # qplot(dat$ADJ_ZSTAT, dat$ADJ_ZSTAT2, alpha=.1) +
    #     ggpubr::stat_cor() +
    #     theme_bw()

    return(dat)
}
