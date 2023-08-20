#' Get the top factors
#'
#' Search a reduction for observations (e.g. traits) that match a given 
#' \code{term} (case-insensitive substring search).
#' Then get the factors with the highest mean loading for the matching features.
#' @param obj \pkg{Seurat} object or dimensionality reduction object.
#' @param obs Phenotype metadata.
#' Not needed if \code{obj} is a \pkg{Seurat} object.
#' @param terms Terms with which to perform substring search of observations.
#' @param search_col Which column in the observation metadata to perform 
#' substring search.
#' @param n_quantiles How many quantiles to bin factor loadings into.
#' @param select_quantiles Which quantiles to return. Defaults to the
#'  top quantile only.
#' @param plot_hist Whether to plot the distribution of loadings.
#' @param verbose Print messages.
#' @inheritParams scKirby::get_obsm
#' @returns data.table
#' 
#' @export
#' @import scKirby
#' @import data.table
#' @importFrom graphics hist
#' @examples
#' obj <- get_HPO()
#' top_factors <- get_top_factors(obj = obj,
#'                                terms = c("parkinson","cardio"),
#'                                search_col = "HPO_label")
get_top_factors <- function(obj,
                            terms,
                            obs = NULL,
                            keys = NULL,
                            search_col = NULL,
                            n_quantiles = 10,
                            select_quantiles = n_quantiles,
                            plot_hist = FALSE,
                            verbose = TRUE) { 
    # devoptera::args2vars(get_top_factors) 
    
    #### Get obsm (embedding) ####
    obsm <- scKirby::get_obsm(obj = obj,
                              keys = keys,
                              verbose = verbose)
    keys <- names(obsm)
    if(length(obsm)>1){
        messager(">1 obsm embedding identified.",
                 "Using first one only:",shQuote(keys[1]),v=verbose)
    }
    obsm <- obsm[[1]]
    #### Get obs (observation metadata) ####
    if (is.null(obs)) obs <- scKirby::get_obs(obj = obj)
    #### Check search_col is present ####
    if (is.null(search_col) | 
        !search_col %in% colnames(obs) | is.null(search_col)) {
        messager("search_col not in obs. Using obs names instead.",v=verbose)
        obs_vec <- scKirby::get_obs_names(obj)
    } else {
        obs_vec <- obs[[search_col]]
    } 
    #### Search obs for matches ####
    select_cols <- rownames(
        obs[grepl(paste(terms,collapse = "|"),obs_vec,ignore.case = TRUE), ]
    )
    #### Report matches ####
    if (length(select_cols) > 0) {
        messager("+", formatC(length(select_cols),big.mark = ","),
                 "matching phenotypes identified:",
                 paste("\n -", select_cols, collapse = ""),v=verbose)
    } else {
        stopper("0 matching phenotypes identified.")
    } 
    #### Handle multiple phenotypes per term
    obsm_means <- if (scKirby::is_class(obsm,"matrix")){
        colMeans(obsm[select_cols,], na.rm = TRUE)
    } else {
        obsm
    }
    if (isTRUE(plot_hist)) {
        graphics::hist(obsm_means, 50, 
                       main = paste("Top", keys[1], "factors:", 
                                    paste(terms,collapse = "; "))) |>
            methods::show()
    }
    #### Filter by quantiles ####
    messager("Filtering by quantiles computed from absolute",
             shQuote(keys[1]),"loadings.",v=verbose)
    quantiles <- cut(obsm_means,
                     breaks = n_quantiles,
                     labels = seq(n_quantiles)) 
    top_factors <- obsm_means[quantiles %in% select_quantiles]
    messager("Returning",formatC(length(top_factors),big.mark = ","),
             "top factor(s).",v=verbose)
    return(top_factors)
}
