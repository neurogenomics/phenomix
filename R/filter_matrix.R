filter_matrix <- function(X,
                          coln_min = 1,
                          rown_min = 1,
                          colsum_min = 0,
                          rowsum_min = 0,
                          features = NULL){ 
    ## Only include subset of features
    if(!is.null(features)){ 
        X <- X[features[features %in% rownames(X)],,drop=FALSE]
    }
    ## Filter out columns with all zeros
    if(!is.null(colsum_min)){
        X <- X[,unname(Matrix::colSums(abs(X), na.rm = TRUE)>=colsum_min),drop=FALSE]  
    } 
    ## Filter out rows with all zeros
    if(!is.null(rowsum_min)){
        X <- X[as.numeric(Matrix::rowSums(abs(X), na.rm = TRUE))>=rowsum_min,,drop=FALSE] 
    }  
    ## Filter columns by number of non-zero entries
    if(!is.null(rown_min)){
        X_nfeat <- apply(X,2, function(x) sum(abs(x)>0, na.rm = TRUE))  
        X <- X[,unname(X_nfeat)>=rown_min,drop=FALSE]  
    }
    ## Filter rows by number of non-zero entries
    if(!is.null(coln_min)){
        X_nsamp <- apply(X,1, function(x) sum(abs(x)>0, na.rm = TRUE))  
        X <- X[unname(X_nsamp)>=coln_min,,drop=FALSE]  
    }
    return(X)
}
