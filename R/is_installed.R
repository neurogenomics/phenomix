is_installed <- function(pkg){
    pkg %in% rownames(installed.packages()) 
}