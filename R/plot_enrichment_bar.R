#' Plot enrichment results: barplot
#'
#' @keywords internal
#' @importFrom dplyr desc arrange
plot_enrichment_bar <- function(res,
                                qvalue_thresh = .05,
                                show_plot = TRUE,
                                save_path = NULL,
                                width = NA,
                                height = NA,
                                dpi = 300,
                                verbose = TRUE) {
    qvalue <- pvalue <- term <- NULL;
    plot_dat <- subset(res, qvalue < .05) |>
        dplyr::arrange(dplyr::desc(pvalue))
    plot_dat$term <- factor(plot_dat$term,
        levels = unique(plot_dat$term),
        ordered = TRUE
    ) 
    gg_bar <- ggplot(
        plot_dat,
        aes(
            y = term, x = -log10(pvalue),
            fill = -log10(pvalue)
        )
    ) +
        geom_bar(stat = "identity") +
        facet_wrap(
            facets = trait ~ .,
            scales = "free_x"
        ) +
        theme_bw()
    if (show_plot) print(gg_bar)
    save_plot(
        save_path = save_path,
        plot = gg_bar,
        width = width,
        height = height,
        dpi = dpi
    )
    return(gg_bar)
}
