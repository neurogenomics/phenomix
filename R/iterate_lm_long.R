iterate_lm_long <- function(xmat,
                            ymat, 
                            BPPARAM, 
                            test_method,
                            multivariate,
                            scale_fn,
                            ...){
    progressbar <- BPPARAM$progressbar
    #### Select lm function ####
    lm_fun <- if(test_method=="glm"){
        iterate_lm_long_glm
    } else if(test_method=="anova") {
        iterate_lm_long_anova
    } else if(test_method=="lm.ridge"){ 
        iterate_lm_long_ridge
    } else if(test_method=="rlm"){ 
        iterate_lm_long_rlm
    } else {
        stopper("test_method not recognized.")
    }
    #### Run ####
    BiocParallel::bplapply(
        BPPARAM = BPPARAM,
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
            res <- lm_fun(dt=dt,
                          multivariate=multivariate,
                          progressbar=progressbar,
                          scale_fn=scale_fn,
                          i=i,
                          ...)
            return(res)
        }) |> 
        data.table::rbindlist(idcol = "yvar",
                              fill = TRUE) 
}