#' Run UMAP
#'
#' Run Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP).
#'
#' Uses \link[uwot]{umap}, but runs and returns PCA by default.
#'
#' @param mat Matrix to run UMAP on.
#' @param transpose Whether to transpose the matrix first.
#' @param add_names Add colnames and rownames to embeddings and loadings.
#' @param seed Seed passed to \[base]{set.seed} for reproducibility between runs.
#' @param ... Additional parameters passed to \link[uwot]{umap}.
#' @inheritParams uwot::umap
#'
#' @importFrom Matrix t
#' @importFrom uwot umap
#' @importFrom dplyr %>%
#'
#' @source \href{https://umap-learn.readthedocs.io/en/latest/}{UMAP documentation}
#' @export
run_umap <- function(mat,
                     transpose = TRUE,
                     pca = ncol(mat),
                     add_names = TRUE,
                     n_components = 2,
                     n_neighbors = 15,
                     min_dist = 0.01,
                     metric = "euclidean",
                     init = "spectral",
                     seed = 2020,
                     verbose = TRUE,
                     ...) {
    if (!is.null(seed)) set.seed(seed)
    #### Check pca arg ####
    if (any(pca > ncol(mat))) {
        pca <- min(pca, ncol(mat), na.rm = TRUE)
        messager(
            "Number of PCA dimensions cannot be",
            "larger than the number of columns in mat.",
            "Setting pca to ", pca
        )
    }
    if (transpose) {
        mat <- Matrix::t(mat)
    }
    umap <- uwot::umap(
        X = as.matrix(mat),
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
    if (add_names) {
        umap$embedding <- umap$embedding %>%
            `colnames<-`(paste0("UMAP.", seq(1, ncol(umap$embedding)))) %>%
            `row.names<-`(rownames(mat))
        for (i in length(umap$pca_models)) {
            pc_embed <- umap$pca_models[[i]]$rotation
            umap$pca_models[[i]]$rotation <- pc_embed %>%
                `colnames<-`(paste0("PC", seq(1, ncol(pc_embed)))) %>%
                `row.names<-`(colnames(mat))
        }
    }
    return(umap)
}
