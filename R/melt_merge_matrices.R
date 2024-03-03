melt_merge_matrices <- function(xmat,
                                ymat){
    ydt <- data.table::as.data.table(ymat,
                                     keep.rownames = "feature") |>
        data.table::setnames(2,"y")
    xdt <- melt_matrix(xmat, 
                       variable.name = "xvar",
                       value.name = "x")
    dt <- xdt[ydt, on="feature"] 
    return(dt)
}