iterate_lm_long <- function(xmat,
                            ymat, 
                            cores, 
                            test_method,
                            multivariate,
                            scale_fn,
                            ...){
    x <- y <- xvar <- NULL; 
    progressbar <- cores$params$progressbar
    add_model_id <- function(res,i){ 
        mid <- gsub("file",paste0("model",i,"_"),basename(tempfile()))
        res[,model_id:=mid,]  
    }
    BiocParallel::bplapply(
        BPPARAM = cores$params,
        X = stats::setNames(seq(ncol(ymat)), 
                            colnames(ymat)), 
        FUN = function(i) {
            tt <- colnames(ymat)[i]
            if(!progressbar){
                messager("-",tt,": (",i,"/",ncol(xmat),")", 
                         parallel=TRUE)    
            }  
            #### Long format ####
            ## Prepare data for rstatix (long format)
            dt <- melt_merge_matrices(xmat = xmat,
                                      ymat = ymat[,tt, drop=FALSE]) 
            dt <- dt[!is.na(x) & !is.na(y),]
            if(isFALSE(dt_var_check(dt,"feature",verbose=!progressbar)) ||
               isFALSE(dt_var_check(dt, "x", verbose=!progressbar)) ||
               isFALSE(dt_var_check(dt, "y", verbose=!progressbar)) ){
                return(NULL)
            } 
            if(test_method=="glm"){
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
                        messager("test_method: glm (univariate)",v=!progressbar)
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
            } else if(test_method=="anova") {
                #### ANOVA ####
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
            } else if(test_method=="lm.ridge"){ 
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
            } else if(test_method=="rlm"){ 
                #### rlm ####
                if(!is.null(scale_fn)){
                    dt[,x:=scale_fn(x)]
                    dt[,y:=scale_fn(y)]
                }
                mod <- MASS::rlm(formula= y~x+xvar,
                                 data=dt,
                                 ...)
                res <- broom::tidy(mod) |>
                    data.table::data.table() 
                add_model_id(res,i)
            } 
            return(res)
        }) |> 
        data.table::rbindlist(idcol = "yvar",
                              fill = TRUE) 
}