#' Plot enrichment results
#'
#' Plot results from \link[phenomix]{iterate_lm} or
#' \link[phenomix]{iterate_gsea}.
#'
#' @param res Enrichment results table.
#' @param types A list of plot types to create
#' @param qvalue_thresh Threshold to filter results.
#' @param top_n Number of top results to plot
#' (sorted by smallest p-values).
#' @param show_plot Print the plot.
#' @param save_path Path to save plot to.
#' @param verbose Print messages.
#' @param ... Additional arguments passed to \link[heatmaply]{heatmaply}.
#'
#' @returns A list of ggplot and/or plotly objects.
#'
#' @export
plot_enrichment <- function(res,
                            types = c("heat", "bar"),
                            qvalue_thresh = .05,
                            top_n = 50,
                            show_plot = TRUE,
                            save_path=NULL,
                            verbose = TRUE,
                            ...) {
    requireNamespace("ggplot2")
    types <- tolower(types)
    res_list <- list()
    if ("bar" %in% types) {
        gg_bar <- plot_enrichment_bar(
            res = res,
            qvalue_thresh = qvalue_thresh,
            show_plot = show_plot,
            save_path = save_path,
            verbose = verbose
        )
        res_list[["barplot"]] <- gg_bar
    }
    if ("heat" %in% types) {
        gg_heat <- plot_enrichment_heat(res,
            qvalue_thresh = .05,
            show_plot = TRUE,
            save_path = NULL,
            verbose = TRUE,
            ...
        )
        res_list[["heatmap"]] <- gg_heat
    }
    return(res_list)
}
