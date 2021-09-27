extract_matrix <- function(obj,
                           assay = NULL,
                           slot = NULL,
                           verbose = TRUE) {
    if (methods::is(obj, "Seurat")) {
        messager("Extracting matrix from Seurat object.", v = verbose)
        mat <- Seurat::GetAssayData(
            object = obj,
            assay = assay,
            slot = slot
        )
        return(mat)
    } else if (is_sparse_matrix(X = obj) |
        methods::is(obj, "matrix") |
        methods::is(obj, "Matrix")) {
        messager("Returning matrix obj directly.", v = verbose)
        return(obj)
    } else {
        stop("obj must be a matrix or Seurat object.")
    }
}
