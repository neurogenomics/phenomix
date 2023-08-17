zscore <- function(dat,
                   column = "BETA",
                   verbose = TRUE) {
    ZSTAT <- NULL;
    column <- column[1]
    
    if(column %in% names(dat)){
        messager("Computing ZSTAT from", paste0(shQuote(column), "."), 
                 v = verbose)
        dat[, ZSTAT:=scale(get(column))]
    } 
}
