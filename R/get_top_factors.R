#' Get the top factors
#'
#' Search a reduction for observations (e.g. traits) that match a given \code{term} (case-insensitive substring search).
#' Then get the factors with the highest mean loading for the matching features.
#'
#' @param obj \pkg{Seurat} object or dimensionality reduction object.
#' @param metadata Phenotype metadata.
#' Not needed if \code{obj} is a \pkg{Seurat} object.
#' @param term Term with which to perform substring search of observations.
#' @param search_col Which column in the observation metadata to perform substring search.
#' @param n_quantiles How many quantiles to bin factor loadings into.
#' @param select_quantiles Which quantiles to return. Defaults to the top quantile only.
#' @param plot_hist Whether to plot the distribution of loadings.
#' @param verbose Print messages.
#' @inheritParams extract_embeddings
#'
#' @return data.table
#' @export
#' @importFrom Seurat Reductions
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% rename group_by slice_max
#' @importFrom data.table data.table
#' @examples
#' data("DEGAS_seurat")
#' top_factors <- get_top_factors(
#'     obj = DEGAS_seurat,
#'     term = "parkinson",
#'     select_quantiles = 8:10
#' )
get_top_factors <- function(obj,
                            metadata = NULL,
                            reduction = NULL,
                            term,
                            search_col = "label_phe",
                            n_quantiles = 10,
                            select_quantiles = n_quantiles,
                            plot_hist = FALSE,
                            verbose = TRUE) {
    embeddings <- extract_embeddings(
        obj = obj,
        reduction = reduction,
        verbose = verbose
    )
    if (is.null(metadata)) metadata <- extract_metadata(obj = obj)
    if (!search_col %in% colnames(metadata) | is.null(search_col)) {
        message(search_col, " not in metadata. Using ", colnames(metadata)[1], " instead.")
        search_col <- colnames(metadata)[1]
    }

    select_cols <- rownames(metadata[grepl(term, unname(metadata[[search_col]]), ignore.case = T), ])
    if (length(select_cols) > 0) {
        message(
            "+ ", length(select_cols), " matching phenotypes identified:\n",
            paste("  -", select_cols, collapse = "\n")
        )
    } else {
        stop("0 matching phenotypes identified.")
    }

    #### Handle multiple phenotypes per term
    if (!is(embeddings, "numeric")) embeddings <- colMeans(embeddings, na.rm = T)
    if (plot_hist) print(hist(embeddings, 50, main = paste("Top", reduction, "factors:", term)))

    quantiles <- cut(abs(embeddings), breaks = n_quantiles, labels = 1:n_quantiles)
    top_factors <- embeddings[quantiles %in% select_quantiles]
    return(top_factors)
}
