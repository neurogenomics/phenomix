rm_empty <- function(res){
    res[!mapply(function(x){is.null(x)||nrow(x)==0},res)]
}