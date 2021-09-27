#' Save ggplot
#'
#' @param save_path Path to save plot to.
#' @param plot ggplot object
#' @param width Plot width
#' @param height Plot height.
#' @param dpi Plot resolution.
#'
#' @keywords internal
#' @importFrom ggplot2 ggsave
save_plot <- function(save_path,
                      plot,
                      width = NA,
                      height = NA,
                      dpi = 300) {
    if (!is.null(save_path)) {
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        ggplot2::ggsave(
            filename = save_path,
            plot = plot,
            width = width,
            height = height,
            dpi = dpi
        )
    }
}
