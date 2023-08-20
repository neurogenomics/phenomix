### From EWCE::generate_celltype_data()
calc_specificity <- function(X) {
    normalised_meanExp <- Matrix::t(Matrix::t(X) * (1 / colSums(X)))
    spec <- normalised_meanExp / (apply(normalised_meanExp, 1, sum) + 1e-12)
    return(spec)
}
