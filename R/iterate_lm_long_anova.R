iterate_lm_long_anova <- function(dt,
                                  progressbar,
                                  scale_fn,
                                  i,
                                  ...){
    x <- y <- xvar <- NULL; 
    if(!is.null(scale_fn)){
        dt[,x:=scale_fn(x)]
        dt[,y:=scale_fn(y)]
    }
    messager("test_method: ANOVA",v=!progressbar) 
    res <- dt |>
        rstatix::group_by(xvar) |>
        rstatix::anova_test(formula = y ~ x,
                            ...) |>
        data.table::data.table()
    add_model_id(res,i)
    return(res)
}