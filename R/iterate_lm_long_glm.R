iterate_lm_long_glm <- function(dt,
                                multivariate,
                                progressbar,
                                scale_fn,
                                i,
                                ...){
    x <- y <- xvar <- NULL; 
    if(isTRUE(multivariate)){
        #### glm: multivariate ####
        messager("test_method: glm (multivariate)",v=!progressbar)
        if(!is.null(scale_fn)){
            dt[,x:=scale_fn(x)]
            dt[,y:=scale_fn(y)]
        }
        mod <- stats::glm(data = dt,
                          formula = y~x*xvar,
                          ...)
        res <- broom::tidy(mod) |>
            data.table::data.table()  
        add_model_id(res,i)
    } else {
        #### glm: univariate ####
        res <- lapply(stats::setNames(unique(dt$xvar),
                                      unique(dt$xvar)),
                      function(xv){
                          messager("test_method: glm (univariate)",
                                   v=!progressbar)
                          dt_sub <- dt[xvar==xv]
                          if(!is.null(scale_fn)){
                              dt_sub[,x:=scale_fn(x)]
                              dt_sub[,y:=scale_fn(y)]
                          }
                          mod <- stats::glm(data = dt_sub,
                                            formula = y~x,
                                            ...) 
                          res <- broom::tidy(mod) |>
                              data.table::data.table() 
                          add_model_id(res,i)
                      }) |> 
            data.table::rbindlist(idcol = "xvar",
                                  fill = TRUE)
    }
    return(res)
}