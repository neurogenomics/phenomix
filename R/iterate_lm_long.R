iterate_lm_long <- function(xmat,
                            ymat, 
                            cores, 
                            method,
                            multivariate,
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
            ## Run tests: glm   
            if(method=="glm"){
                if(isTRUE(multivariate)){
                    messager("Method: glm (multivariate)",v=!progressbar)
                    mod <- stats::glm(data = dt,
                                      formula = y~x*xvar,
                                      ...)
                    res <- broom::tidy(mod) |>
                        data.table::data.table()  
                    add_model_id(res,i)
                } else {
                    res <- lapply(stats::setNames(unique(dt$xvar),
                                                  unique(dt$xvar)),
                                  function(xv){
                        messager("Method: glm (univariate)",v=!progressbar)
                        mod <- stats::glm(data = dt[xvar==xv],
                                          formula = y~x,
                                          ...)
                        res <- broom::tidy(mod) |>
                            data.table::data.table() 
                        add_model_id(res,i)
                    }) |> 
                        data.table::rbindlist(idcol = "xvar",
                                              fill = TRUE)
                }
            } else if(method=="anova") {
            ## Run tests: ANOVA 
                messager("Method: ANOVA",v=!progressbar) 
                res <- dt |>
                rstatix::group_by(xvar) |>
                rstatix::anova_test(formula = y ~ x,
                                    ...) |>
                data.table::data.table()
            } else if(method=="lm.ridge"){ 
            ### XGboost ####
                mod <- MASS::lm.ridge(formula= y~x+xvar,
                                      data=dt,
                                      ...)
                res <- broom::tidy(mod) |>
                    data.table::data.table() 
                add_model_id(res,i)
            } else if(method=="rlm"){ 
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