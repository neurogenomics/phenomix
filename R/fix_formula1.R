fix_formula1 <- function(formula,
                         dat){
    
    formula_split <- as.character(formula)
    predictors <- trimws(strsplit(tail(formula_split,1),"[+]")[[1]])
    predictors1 <- predictors[predictors %in% colnames(dat)]
    formula1 <- formula(paste(formula_split[2],
                              formula_split[1],
                              paste(predictors1,collapse = " + "))) 
    return(formula1)
}