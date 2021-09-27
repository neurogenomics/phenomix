#' Plot top trait-trait correlations
#'
#' Takes the output of \link[phenomix]{find_neighbors} as input to \code{knn}.
#'
#' @param knn A melted similarity matrix or K-nearest neighbor graph.
#' @param top_n The number of top correlations to plot.
#' @param non_self Remove trait2 that are present in trait1.
#' @param show_plot Whethe to print the plot.
#'
#' @return \code{ggplot} object.
#'
#' @export
#' @import ggplot2
plot_trait_cor <- function(knn,
                           top_n = 10,
                           non_self = TRUE,
                           show_plot = TRUE) {
    if (non_self) {
        plot_dat <- subset(knn, !trait2 %in% trait1)
    } else {
        plot_dat <- knn
    }
    if (!is.null(top_n)) {
        plot_dat <- plot_dat[seq(1, top_n)]
    }

    gg_bar <- ggplot(
        plot_dat,
        aes(x = trait1, y = similarity, fill = similarity)
    ) +
        geom_bar(stat = "identity") +
        facet_grid(facets = trait2 ~ .) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "white"),
            strip.text.y = element_text(angle = 0)
        )
    if (show_plot) print(gg_bar)
    return(gg_bar)
}
