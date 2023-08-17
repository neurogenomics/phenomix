check_required_cols <- function(dat,
                                 cols=c("CHR","BP","SNP")){ 
    missing <- cols[!cols %in% names(dat)]
    if(length(missing)>0){
        stopper("Missing required columns:",
                paste("\n -",missing, collapse = ""))
    }
}