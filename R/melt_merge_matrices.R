melt_merge_matrices <- function(xmat,
                                ymat){
    xdt <- data.table::as.data.table(xmat,
                                     keep.rownames = "feature") |>
        data.table::setnames(2,"x")
    ydt <- melt_matrix(ymat, 
                       variable.name = "yvar",
                       value.name = "y")
    dt <- xdt[ydt, on="feature"] 
    return(dt)
}