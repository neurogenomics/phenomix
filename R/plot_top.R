#' Plot the top factors
#'
#' Plot the top factors for a given set of phenotypes,
#' OR plot the top features for a given set of factors.
#'
#' @param top_data Data with the top phenotypes or features.
#' @param factors Which factors to plot.
#' @param x x-axis variable.
#' @param y y-axis variable.
#' @param fill fill variable.
#' @param title Plot title.
#' @param invert_vars Switch the axes of factors and loadings in plot.
#' @param show_plot Whether to print the plot or simply return it.
#'
#' @return ggplot object
#' @keywords internal 
#' @importFrom stringr str_split
#' @importFrom methods is
plot_top <- function(top_data,
                     factors = NULL,
                     x = "feature",
                     y = "loading",
                     fill = "factor",
                     title = NULL,
                     invert_vars = FALSE,
                     show_plot = TRUE) {
    # devoptera::args2vars(plot_top)
    requireNamespace("ggplot2")
    
    if (is.null(factors)) factors <- seq(length(unique(top_data$factor)))
    if (methods::is(factors, "character")) {
        factors <- stringr::str_split(factors, "_", simplify = TRUE)[, 2]|>
            as.numeric()
    }
    if (isTRUE(invert_vars)) {
        x1 <- x
        fill1 <- fill
        x <- fill
        fill <- x1
    }
    factor_names <- unique(top_data$factor)[factors]
    #### Filtering  ####
    top_data <- top_data[top_data[["factor"]] %in% factor_names, ]
    #### Plot ####
    gp <- ggplot(top_data, 
                 mapping = aes_string(x = x, y = y, fill = fill)) +
        geom_bar(stat = "identity") +
        coord_polar() +
        theme_minimal() +
        labs(title = title)
    if (isTRUE(show_plot)) methods::show(gp)
    return(gp)
}
