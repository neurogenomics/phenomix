normalise <- function(mat,
                      method = "log1p",
                      verbose = TRUE,
                      ...){
    opts <- c("log1p","yeojohnson")
    method <- tolower(method)[1]
    if(method=="log1p"){
        messager("Normalising method: log1p",v=verbose)
        mat_norm <- log1p(mat)
    } else if (method %in% c("yj","yeojohnson")){
        messager("Normalising method: YeoJohnson",v=verbose)
        mat_norm <- normalise_yeojohnson(mat = mat,
                                         prop = 3/4,
                                         ...) 
    } else {
        stop_msg <- paste0("Normalisation must be one of:\n",
                          paste0(" - ",opts, collapse = "\n"))
        stop(stop_msg)
    }
    return(mat_norm)
}