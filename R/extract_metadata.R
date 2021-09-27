extract_metadata <- function(obj) {
    # obj <- DEGAS_contributionGene; obj <- scNLP::pseudo_seurat
    if (is(obj, "DimReduc")) {
        message("Using embedding rownames as metadata")
        metadata <- data.frame(
            id = rownames(obj@cell.embeddings),
            ### Redundant but extra column prevents df from turning into list sometimes
            label_phe_code = rownames(obj@cell.embeddings),
            row.names = rownames(obj@cell.embeddings)
        )
        return(metadata)
    }
    if (is(obj, "Seurat")) {
        metadata <- obj@meta.data
        return(metadata)
    }
}
