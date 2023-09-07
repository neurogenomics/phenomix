#' Plot the top feature loadings
#'
#' Plot the top feature loadings for from a given reduction. 
#' @param obj \pkg{Seurat} object.
#' @param n_features Number of top features per factor to select.
#' @param n_features_plot Number of top features to plot. 
#' @inheritParams plot_top
#' @inheritParams scKirby::get_varm
#'
#' @return top_features \code{data.table}
#' @export
#' @import scKirby
#' @import data.table
#' @examples
#' obj <- get_HPO()
#' top_features <- get_top_features(obj = obj)
get_top_features <- function(obj,
                             keys = NULL,
                             n_features = 3,
                             n_features_plot = 20,
                             factors = seq(10),
                             invert_vars = FALSE,
                             fill = "factor",
                             show_plot = TRUE,
                             title = NULL,
                             verbose = TRUE) {
    # devoptera::args2vars(get_top_features)
    
    loading <- NULL;
    varm <- scKirby::get_varm(obj = obj,
                              keys = keys,
                              verbose = verbose)
    keys <- names(varm)
    varm <- get_one_element(l = varm, 
                            verbose = verbose) 
    top_features <- (
        data.table::as.data.table(varm,keep.rownames = "feature") |>
        data.table::melt.data.table(id.vars = "feature", 
                                    variable.name = "factor", 
                                    value.name = "loading")
    )[,.SD[abs(loading) %in% head(sort(abs(loading)), n_features)],by="factor"]
     
    if (isTRUE(show_plot)) {
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
