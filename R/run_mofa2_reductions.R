run_mofa2_reductions <- function(model,
                                 reductions,
                                 verbose = TRUE) {
    for (r in tolower(gsub("-|_", "", reductions))) {
        if (r == "umap") {
            messager("Running UMAP on MOFA2 factors.", v = verbose)
            try({
                model <- MOFA2::run_umap(model)
            })
        }
        if (r == "tsne") {
            messager("Running t-SNE on MOFA2 factors.", v = verbose)
            try({
                model <- MOFA2::run_tsne(model)
            })
        }
    }
    return(model)
}
