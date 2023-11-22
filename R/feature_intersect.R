feature_intersect <- function(X_list,
                              verbose=TRUE){
    rn <- lapply(X_list,rownames)
    f <- Reduce(intersect,rn)
    if(length(f)==0) {
        stopper("No intersecting features between",length(X_list),
                "matrices.",v=verbose)
    }
    messager(format(length(f), big.mark = ","), 
            "intersecting features between all matrices.",v=verbose)
    return(f)
}