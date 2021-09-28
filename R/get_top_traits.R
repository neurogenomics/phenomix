#' Plot the top traits
#'
#' Create a radar chart of the traits (phenotypes) with the
#' highest loadings for each reduction factor.
#'
#' @param obj \pkg{Seurat} object or dimensionality reduction object.
#' @param metadata Phenotype metadata.
#' Not needed if \code{obj} is a \pkg{Seurat} object.
#' @param n_traits Number of top traits per reduction factor to select.
#' @param verbose Print messages.
#' @inheritParams plot_top
#' @inheritParams extract_embeddings
#'
#' @return \code{data.table} of top traits.
#' 
#' @export
#' @importFrom Seurat Reductions
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% rename group_by slice_max
#' @importFrom data.table data.table
#' @examples
#' degas <- get_DEGAS()
#' top_phenos <- get_top_traits(obj = degas)
get_top_traits <- function(obj,
                           metadata = NULL,
                           reduction = NULL,
                           n_traits = 3,
                           factors = seq(1,10),
                           invert_vars = FALSE,
                           show_plot = TRUE,
                           title = NULL,
                           x = "phenotype",
                           y = "loading",
                           verbose = TRUE) {
    Var1 <- Var2 <- value <- loading <- NULL;
    embeddings <- extract_embeddings(
        obj = obj,
        reduction = reduction,
        verbose = verbose
    )
    if (is.null(metadata)) metadata <- extract_metadata(obj = obj)
    top_traits <- reshape2::melt(embeddings) %>%
        merge(metadata,
            by.x = "Var1", by.y = 0, all.x = TRUE
        ) %>%
        dplyr::rename(phenotype = Var1,
                      factor = Var2, 
                      loading = value) %>%
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading),
                         n = n_traits) %>%
        data.table::data.table()

    if (show_plot) {
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
