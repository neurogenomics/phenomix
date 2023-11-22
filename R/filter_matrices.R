filter_matrices <- function(X_list,
                            coln_min = 1,
                            rown_min = 1,
                            colsum_min = 0,
                            rowsum_min = 0,
                            verbose = TRUE){
    
    features <- feature_intersect(X_list,
                                  verbose = verbose) 
    lapply(X_list, function(x){
        filter_matrix(x, 
                      coln_min = coln_min,
                      rown_min = rown_min,
                      colsum_min = colsum_min,
                      rowsum_min = rowsum_min,
                      features = features)
    })
}