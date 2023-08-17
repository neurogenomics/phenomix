add_mofa2_dimred <- function(obj,
                             model,
                             assay = NULL,
                             reduction = c("mofa2", "umap")[1],
                             key = reduction) {
    reduction <- tolower(reduction[1])
    if (reduction == "umap") {
        key <- paste0("mofa2", key)
        cnames <- model@samples_metadata$sample
        embeddings <- data.frame(
            UMAP1 = model@dim_red$UMAP$UMAP1,
            UMAP2 = model@dim_red$UMAP$UMAP2,
            row.names = scKirby::get_obs_names(obj)
        )
        DR <- Seurat::CreateDimReducObject(
            embeddings = as.matrix(embeddings),
            key = paste0(key, "_"),
            assay = assay
        )
        DR@misc[["variance_explained"]] <- model@cache$variance_explained
        obj[[key]] <- DR
    }
    if (reduction == "mofa2") {
        cnames <- model@samples_metadata$sample
        embeddings <- model@expectations$Z[[1]]
        # MOFA2 modifies the names. Restore them to original names here
        rownames(embeddings) <- colnames(obj)
        loadings <- model@expectations$W[[1]]
        # model@cache$variance_explained
        DR <- Seurat::CreateDimReducObject(
            embeddings = as.matrix(embeddings),
            loadings = as.matrix(loadings),
            key = paste0(key, "_"),
            assay = assay
        )
        DR@misc[["variance_explained"]] <- model@cache$variance_explained
        obj[[key]] <- DR
    }

    return(obj)
}
