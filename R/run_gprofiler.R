#' Run gprofiler
#'
#' Run \link[gprofiler2]{gost} to perform
#'  Gene Set  Enrichment Analysis (GSEA).
#'
#' @param factors Which factors to run gprofiler on.
#' @param verbose Print messages.
#' @param ... Additional arguments passed to \link[gprofiler2]{gost}.
#' @inheritParams get_top_features
#' @inheritParams gprofiler2::gost
#'
#' @export
#' @importFrom gprofiler2 gost
run_gprofiler <- function(obj,
                          n_features = 3,
                          factors = seq(1, 4),
                          keys = NULL,
                          show_plot = TRUE,
                          verbose = TRUE,
                          ...) {
    top_genes <- get_top_features(
        obj = obj,
        keys = keys,
        n_features = n_features,
        show_plot = FALSE,
        verbose = verbose
    )$data
    messager("Running gprofiler on", length(factors), "factors",
        v = verbose
    )
    if (!is.character(factors)) {
        factors <- as.character(sort(unique(top_genes$factor))[factors])
    }

    gene_lists <- lapply(factors, function(f) {
        as.character(subset(top_genes, factor %in% f)$feature)
    }) |> `names<-`(factors)

    gres <- gprofiler2::gost(query = gene_lists, ...)
    if (isTRUE(show_plot)) {
        gg_gprof <- gprofiler2::gostplot(gres)
        print(gg_gprof)
    }
    res <- data.table::data.table(gres$result)
    return(res)
}
