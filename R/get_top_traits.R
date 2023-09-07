#' Plot the top traits
#'
#' Create a radar chart of the traits (phenotypes) with the
#' highest loadings for each reduction factor.
#'
#' @param obj \pkg{Seurat} object or dimensionality reduction object.
#' @param obs Phenotype metadata.
#' Not needed if \code{obj} is a \pkg{Seurat} object.
#' @param n_traits Number of top traits per reduction factor to select.
#' @param verbose Print messages.
#' @inheritParams plot_top
#' @inheritParams scKirby::get_obsm
#' @returns \link[data.table]{data.table} of top traits.
#' 
#' @export 
#' @import data.table
#' @examples
#' obj <- get_HPO()
#' top_phenos <- get_top_traits(obj = obj)
get_top_traits <- function(obj,
                           obs = NULL,
                           keys = NULL,
                           n_traits = 3,
                           factors = seq(10),
                           invert_vars = FALSE,
                           show_plot = TRUE,
                           title = NULL,
                           x = "trait",
                           y = "loading",
                           verbose = TRUE) {
    # devoptera::args2vars(get_top_traits)
    
    loading <- NULL;
    obsm <- scKirby::get_obsm(obj = obj,
                              keys = keys,
                              verbose = verbose)
    keys <- names(obsm)
    obsm <- get_one_element(l = obsm, 
                            verbose = verbose) 
    if (is.null(obs)) obs <- scKirby::get_obs(obj = obj)
    
    top_traits <- (
        data.table::as.data.table(obsm, keep.rownames = "trait") |>
            data.table::melt.data.table(id.vars = "trait", 
                                        variable.name = "factor", 
                                        value.name = "loading")
    )[,.SD[abs(loading) %in% head(sort(abs(loading)), n_traits)],
      by="factor"] |>
        merge(data.table::as.data.table(obs,keep.rownames = "trait"),
              by="trait") 
    #### Plot ####
    if (isTRUE(show_plot)) {
        gp <- plot_top(
            top_data = top_traits,
            x = x,
            y = y,
            factors = factors,
            invert_vars = invert_vars,
            title = title
        )
    } else { gp <- NULL}
    return(list(data=top_traits,
                plot=gp))
}
