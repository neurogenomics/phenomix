map_xref <- function(dat,
                     prefix="MONDO",
                     new_col=paste0(tolower(prefix),"_id"),
                     verbose=TRUE){
    dbXRefs <- id <- NULL
    messager("Adding xref column:",new_col,v=verbose)
    dat[grepl(paste0("^",prefix),id),(new_col):=id]
    dat[get(new_col)=="NA",(new_col):=NA]
    if(sum(is.na(unlist(dat[[new_col]])))>0) {
        dat[is.na(get(new_col)),(new_col):=sapply(dbXRefs,function(x){
            r <- grep(paste0("^",prefix),unlist(x),value = TRUE)
            if(length(r)==0) NA else unlist(r)
        })]
    }
    dat[get(new_col)=="NA",(new_col):=NA]
    messager("% of non-NA rows:",
             round(sum(!is.na(dat[[new_col]]))/nrow(dat)*100,2),v=verbose)
}