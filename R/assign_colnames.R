assign_colnames <- function(obj,
                            cnames) {
    if (methods::is(obj, "Seurat")) {
        obj <- Seurat::RenameCells(
            object = obj,
            new.names = cnames
        )
    } else {
        colnames(obj) <- cnames
    }
    return(obj)
}
