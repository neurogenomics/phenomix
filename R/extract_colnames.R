extract_colnames <- function(obj) {
    if (methods::is(obj, "Seurat")) {
        cnames <- colnames(obj)
    }
    if (methods::is(obj, "list")) {
        cnames <- lapply(seq(1, length(obj)), function(i) {
            colnames(obj[[i]])
        })
    }
    if (methods::is(obj, "data.frame")) {
        cnames <- colnames(obj)
    }
    if (methods::is(obj, "matrix") | methods::is(obj, "Matrix")) {
        cnames <- colnames(obj)
    }
    return(cnames)
}
