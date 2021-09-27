#' Plot the top phenotypes
#'
#' Create a radar chart of the phenotypes with the
#' highest loadings for each reduction factor.
#'
#' @param seurat \code{Seurat} object or \code{Seurat::DimReducObject}.
#' @param reduction Reduction to use. If \code{NULL}, defaults to first available reduction.
#' @param n_phenotypes Number of top phenotypes per reduction factor to select.
#' @inheritParams plot_top
#'
#' @return \code{data.table} of top phenotypes
#' @export
#' @importFrom Seurat Reductions
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% rename group_by slice_max
#' @importFrom data.table data.table
#' @examples
#' data("DEGAS_seurat")
#' top_phenos <- get_top_phenotypes(seurat = DEGAS_seurat)
get_top_phenotypes <- function(seurat,
                               reduction = NULL,
                               n_phenotypes = 3,
                               factors_plot = 1:10,
                               invert_vars = FALSE,
                               show_plot = TRUE,
                               title = NULL,
                               x = "phenotype",
                               y = "loading") {
    loadings <- extract_loadings(obj = seurat, reduction = reduction)
    metadata <- extract_metadata(obj = seurat)
    top_phenos <- reshape2::melt(loadings) %>%
        merge(metadata,
            by.x = "Var1", by.y = 0, all.x = TRUE
        ) %>%
        dplyr::rename(phenotype = Var1, factor = Var2, loading = value) %>%
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading), n = n_phenotypes) %>%
        data.table::data.table()

    if (show_plot) {
        gp <- plot_top(
            top_data = top_phenos,
            x = x,
            y = y,
            factors = factors_plot,
            invert_vars = invert_vars,
            title = title
        )
    }
    return(top_phenos)
}
