#' Iterate GSEA
#' 
#' Iterate gene set enrichment analysis across all traits and celltypes.
#'
#' @param xmat gene x trait matrix.
#' @param ymat gene x celltype matrix.
#' @param correction_method Multiple-testing correction method 
#' to be passed to \code{stats::p.adjust}.
#' @param qvalue_thresh q.value threshold to use when report 
#' significant results summary.
#' @param x_quantiles The number of quantiles to bin \code{ymat} data into.
#' @param y_quantiles The number of quantiles to bin \code{ymat} data into.
#' @param use_quantiles Which quantiles in to use in 
#' \code{xmat} and \code{ymat}.
#' @param nCores Number of cores to use in parallel. 
#' Will optimize if \code{NULL}.
#' 
#' @return \code{data.table} of enrichment results. 
#' 
#' @export
#' @importFrom parallel mclapply
#' @importFrom GeneOverlap newGeneOverlap testGeneOverlap
#' @importFrom broom tidy
#' @importFrom data.table rbindlist
#' @importFrom dplyr mutate
#' @importFrom stats p.adjust
#' @importFrom methods slot
#' @examples
#' ### DeGAs loadings
#' obj <- get_DEGAS()
#' xmat <- get_varm(obj = obj)
#' xmat <- xmat[, 1:10] # Let's use just 10 components as an example
#'
#' ### Celltype Dataset
#' ctd <- get_BlueLake2018_FrontalCortexOnly()
#' ymat <- ctd[[1]]$specificity
#' res_gsea <- iterate_gsea(xmat = xmat,
#'                          ymat = ymat,
#'                          nCores = 1)
iterate_gsea <- function(xmat,
                         ymat,
                         correction_method = "BH",
                         qvalue_thresh = .05,
                         x_quantiles = 10,
                         y_quantiles = 10,
                         use_quantiles = 10,
                         nCores = 1) {
    trait <- celltype <- p.value <- qvalue <- NULL;
    gene_intersect <- intersect(rownames(xmat), rownames(ymat))
    message(length(gene_intersect), 
            " intersecting genes between GWAS and CTD matrices.")
    ### Run lm  for all celltypes against this trait
    messager(
        "Running ", formatC(ncol(xmat) * ncol(ymat), big.mark = ","),
        " tests: ", formatC(ncol(xmat), big.mark = ","), " traits x ",
        formatC(ncol(ymat), big.mark = ","), " celltypes."
    )

    gsea_res <- parallel::mclapply(1:ncol(xmat), function(i) {
        tt <- colnames(xmat)[i]
        messager(" - ", tt, ": (", i, "/", ncol(xmat), ")", parallel=TRUE)
        lapply(colnames(ymat), function(ct) { 
            dat <- data.frame(
                trait = cut(xmat[gene_intersect, tt], 
                            breaks = x_quantiles, 
                            labels = 1:x_quantiles),
                celltype = cut(ymat[gene_intersect, ct], 
                               breaks = y_quantiles, 
                               labels = 1:y_quantiles),
                row.names = gene_intersect
            )
            res <- GeneOverlap::newGeneOverlap(
                listA = rownames(subset(dat, trait %in% use_quantiles)),
                listB = rownames(subset(dat, celltype %in% use_quantiles))
            ) |>
                GeneOverlap::testGeneOverlap()
            res_df <- data.frame(
                term = ct,
                trait_genes = length(res@listA),
                celltype_genes = length(res@listB),
                intersection = length(res@intersection),
                union = length(res@union),
                genome.size = length(res@genome.size),
                odds.ratio = res@odds.ratio,
                Jaccard = res@Jaccard,
                p.value = res@pval
            )
            return(res_df)
        }) |> data.table::rbindlist()
    }, mc.cores = nCores) |>
        `names<-`(colnames(xmat)) |>
        data.table::rbindlist(idcol = "trait")

    ### Multiple-testing correction
    gsea_res <- gsea_res |>
        dplyr::mutate(qvalue = stats::p.adjust(p = p.value,
                                               method = correction_method)) |>
        dplyr::rename(pvalue = p.value)
    ### Filter only sig results
    sig_res <- gsea_res |>
        subset(qvalue < qvalue_thresh)
    messager("\n", formatC(nrow(sig_res), big.mark = ","),
            " significant results @ ", 
            correction_method, "<", qvalue_thresh)
    ### Return FULL results (not just sig)
    return(gsea_res)
}
