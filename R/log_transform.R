log_transform <- function(dat,
                          log_vars,
                          verbose=TRUE){
    for(x in log_vars){
        if(x %in% colnames(dat)){
            messager("Performing natural log-tranformation on",paste0(x,"."),
                     v=verbose)
            dat[,paste0("log",x) :=log(get(x))]
        } 
    } 
}