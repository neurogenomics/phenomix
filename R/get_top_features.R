#' Plot the top feature loadings
#'
#' Plot the top feature loadings for from a given reduction.
#'
#' @param obj \pkg{Seurat} object.
#' @param reduction Reduction to use. If \code{NULL}, defaults to first available reduction.
#' @param n_features Number of top features per factor to select.
#' @param n_features_plot Number of top features to plot.
#' @param verbose Print messages.
#' @inheritParams plot_top
#'
#' @return top_features \code{data.table}
#' @export
#' @importFrom Seurat Reductions
#' @importFrom reshape2 melt
#' @importFrom dplyr rename group_by slice_max
#' @importFrom data.table data.table
#' @examples
#' degas <- get_DEGAS()
#' top_features <- get_top_features(obj = degas)
get_top_features <- function(obj,
                             reduction = NULL,
                             n_features = 3,
                             n_features_plot = 20,
                             factors = seq(1, 10),
                             invert_vars = FALSE,
                             fill = "factor",
                             show_plot = TRUE,
                             title = NULL,
                             verbose = TRUE) {
    Var1 <- Var2 <- value <- loading <- NULL;
    loadings <- get_varm(
        obj = obj,
        reduction = reduction,
        verbose = verbose
    )
    top_features <- reshape2::melt(loadings) |>
        dplyr::rename(feature = Var1, factor = Var2, loading = value) |>
        dplyr::group_by(factor) |>
        dplyr::slice_max(order_by = abs(loading), n = n_features) |>
        data.table::data.table()

    if (show_plot) {
        gp <- plot_top(
            top_data = top_features,
            factors = factors,
            invert_vars = invert_vars,
            x = "feature",
            y = "loading",
            fill = fill,
            title = title
        )
    } else {
        gp <- NULL
    }
    return(list(
        data = top_features,
        plot = gp
    ))
}
