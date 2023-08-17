#' Plot enrichment results: heatmap
#'
#' @keywords internal
#' @importFrom dplyr desc arrange
#' @importFrom data.table data.table dcast.data.table
#' @importFrom tibble column_to_rownames
#' @importFrom heatmaply heatmaply
plot_enrichment_heat <- function(res,
                                 qvalue_thresh = .05,
                                 show_plot = TRUE,
                                 save_path = NULL,
                                 verbose = TRUE,
                                 ...) {
    pvalue <- qvalue <- 
    plot_dat <- subset(res, qvalue < .05) |>
        dplyr::arrange(dplyr::desc(pvalue))
    plot_dat$term <- factor(plot_dat$term,
        levels = unique(plot_dat$term),
        ordered = TRUE
    )

    dmat <- data.table::data.table(plot_dat) |>
        data.table::dcast.data.table(
            formula = term ~ trait,
            value.var = "negLogP",
            fun.aggregate = mean,
            fill = 0
        ) |>
        tibble::column_to_rownames("term") |>
        as.matrix()
    # row_side_colors <- data.frame(stringr::str_split(rownames(dmat),"[.]",
    #                                                  simplify = TRUE)) |>
    #     `colnames<-`(c("Tissue","Celltype"))
    gg_heat <- heatmaply::heatmaply(
        x = dmat,
        ...
    )

    # gplots::heatmap.2(dmat)
    if (show_plot) print(gg_heat)
    return(gg_heat)
}
