iterate_lm_long_rlm <- function(dt,
                                progressbar,
                                scale_fn,
                                i,
                                ...){
    x <- y <- NULL; 
    if(!is.null(scale_fn)){
        dt[,x:=scale_fn(x)]
        dt[,y:=scale_fn(y)]
    }
    messager("test_method: rlm",v=!progressbar) 
    mod <- MASS::rlm(formula= y~x+xvar,
                     data=dt,
                     ...)
    res <- broom::tidy(mod) |>
        data.table::data.table() 
    return(res)
}