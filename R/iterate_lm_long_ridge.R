iterate_lm_long_ridge <- function(dt,
                                  progressbar,
                                  scale_fn,
                                  i,
                                  ...){
    x <- y <- NULL; 
    messager("test_method: lm.ridge",v=!progressbar)
    #### lm.ridge ####
    if(!is.null(scale_fn)){
        dt[,x:=scale_fn(x)]
        dt[,y:=scale_fn(y)]
    }
    mod <- MASS::lm.ridge(formula= y~x+xvar,
                          data=dt,
                          ...)
    res <- broom::tidy(mod) |>
        data.table::data.table() 
    add_model_id(res,i)
    return(res)
}