dt_var_check <- function(dt,
                      var,
                      min_n=2,
                      verbose=TRUE){
    is_ok <- !length(unique(dt[[var]])) < min_n
    if(isFALSE(is_ok)){
        messager("Not enough unique values:",shQuote(var),
                 "Returning NULL.",
                 v=verbose) 
    }
    return(is_ok)
}