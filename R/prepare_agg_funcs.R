prepare_agg_funcs <- function(merged_hits,
                              extra_cols=c("ABC.Score","CellType"),
                              verbose=TRUE){
    #### Define functions ####
    mean_ <- function(x){mean(as.numeric(x),na.rm=T)}
    sd_ <- function(x){sd(as.numeric(x),na.rm=T)}
    uniqueN_ <- function(x){data.table::uniqueN(x, na.rm=TRUE)}
    N_ <- function(x){length(x)}
    collapse_ <- function(x){paste(unique(x),collapse = ";")}
    #### Define columns ####
    agg_functions = list(BETA_mean=mean_, P_mean=mean_,SE_mean=mean_, 
                         BETA_sd=sd_, P_sd=sd_, SE_sd=sd_,
                         nonuniqueSNPS=N_, NSNPS=uniqueN_) 
    col_selection <- c("BETA","P","SE",
                       "BETA","P","SE",
                       "SNP","SNP")
    #### Add extra columns ####
    extra_cols <- extra_cols[extra_cols %in% colnames(merged_hits)]
    if(length(extra_cols)>0){
        messager("Extra cols used:",paste(extra_cols,collapse = ", "),
                 v=verbose)
        list2 <- list()
        if("ABC.Score" %in% extra_cols){
            agg_functions[[paste0("ABC.Score","_mean")]] <- mean_
            col_selection <- c(col_selection, "ABC.Score")
            agg_functions[[paste0("ABC.Score","_sd")]] <- sd_
            col_selection <- c(col_selection, "ABC.Score")
        }
        if("CellType" %in% extra_cols){
            agg_functions[["CellTypes"]] <- collapse_
            col_selection <- c(col_selection, "CellType")
        } 
    }
    return(list(agg_functions=agg_functions,
                col_selection=col_selection))
}
