#' Regress out the effect of gene attributes
#'
#' Queries \code{biomaRt} for \code{gene_attributes} and then regresses them out of \code{xmat}.
#'
#' @param xmat gene x sample matrix.
#' @param attributes Gene attributes to extract from 
#' \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#'  and then regress from \code{xmat}.
#' @inheritParams iterate_lm
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' ctd <- phenomix::get_BlueLake2018_FrontalCortexOnly()
#' xmat <- ctd[[1]]$mean_exp[1:100, ]
#' adjusted_df <- phenomix:::regress_gene_info(xmat)
#' } 
#' @export
regress_gene_info <- function(xmat,
                              attributes = c("GENELEN"),
                              correction_method = "BH",
                              verbose=TRUE) {
    gene_info <- get_gene_length(
        gene_hits = data.table::data.table(gene=rownames(xmat)), 
        gene_var = "gene",
        use_symbols = TRUE,
        verbose = verbose
    )
    gene_intersect <- intersect(rownames(xmat), gene_info$gene)
    xdat <- as.matrix(xmat[gene_intersect, ])
    ydat <- as.matrix(gene_info[gene %in% gene_intersect, get(attributes)])
    #### Run model
    messager("Training model.",v=verbose)
    mod <- stats::lm(xdat ~ ydat)
    messager("Generating predictions.",v=verbose)
    adjusted_df <- stats::predict(mod, data.frame(xdat))
    return(adjusted_df)
}
