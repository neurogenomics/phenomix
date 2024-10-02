#' Run UMAP
#'
#' Run Uniform Manifold Approximation and Projection for Dimension Reduction 
#' (UMAP).
#'
#' Uses \link[uwot]{umap}, but runs and returns PCA by default. 
#' @param seed Seed passed to \link[base]{set.seed} for reproducibility 
#' between runs.
#' @param add_names Add observation (sample) names to the rows of the embedding.
#' @inheritParams scKirby::get_x
#' @inheritParams uwot::umap
#' @inheritDotParams uwot::umap
#'
#' @importFrom Matrix t
#' @importFrom uwot umap 
#'
#' @source \href{https://umap-learn.readthedocs.io/en/latest/}{
#' UMAP documentation}
#' @export
#' @examples
#' obj <- get_HPO()[seq(50),seq(100)]
#' um <- run_umap(obj)
run_umap <- function(obj,
                     transpose = TRUE,
                     pca = NULL,
                     add_names = TRUE,
                     n_components = 2,
                     n_neighbors = 30,
                     min_dist = 0.01,
                     metric = "euclidean",
                     init = "spectral",
                     seed = 2020,
                     verbose = TRUE,
                     ...) {
    # devoptera::args2vars(run_umap)
    
    if (!is.null(seed)) set.seed(seed)
    X <- obj
    # X <- scKirby::get_x(obj = obj,
    #                     n = 1,
    #                     as_sparse = FALSE,
    #                     transpose = transpose,
    #                     verbose = verbose)
    #### Check pca arg ####
    # if(is.null(pca)) pca <- ncol(X) 
    if (any(pca > ncol(X))) {
        pca <- min(pca, ncol(X), na.rm = TRUE)
        messager(
            "Number of PCA dimensions cannot be",
            "larger than the number of columns in X.",
            "Setting pca to", paste0(pca,".")
        )
    }
    #### Check n_neighbors arg ####
    if(n_neighbors>=nrow(X)){
        n_neighbors <- nrow(X)-1
        messager("n_neighbors must be smaller than the dataset size.",
                 paste0("Setting n_neighbors=",n_neighbors),v=verbose)
    }
    #### Run UMAP ####
    um <- uwot::umap(
        X = as.matrix(X), ## Must be a dense matrix in uwot implementation
        ret_model = TRUE,
        ret_nn = TRUE,
        pca = pca,
        n_neighbors = n_neighbors,
        n_components = n_components,
        min_dist = min_dist,
        metric = metric,
        init = init,
        verbose = verbose,
        ...
    ) 
    if (isTRUE(add_names)) {
        um$embedding <- um$embedding |>
            `colnames<-`(paste0("UMAP.", seq(ncol(um$embedding)))) |>
            `row.names<-`(rownames(X))
        tryCatch({
            for (i in length(um$pca_models)) {
                pc_embed <- um$pca_models[[i]]$rotation
                um$pca_models[[i]]$rotation <- pc_embed |>
                    `colnames<-`(paste0("PC", seq(ncol(pc_embed)))) |>
                    `row.names<-`(colnames(X))
            }
        }, error= function(e) message(e))
    }
    return(um)
}
