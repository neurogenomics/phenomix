zscore <- function(dat,
                   column = "BETA",
                   verbose = TRUE) {
    messager("Computing ZSTAT from", paste0(column, "."), v = verbose)
    data.table::setDT(dat)[, ZSTAT := scale(get(column))]
}
