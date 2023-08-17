fix_formula1 <- function(formula,
                         dat,
                         verbose=TRUE) {
    
    formula_split <- as.character(formula)
    predictors <- trimws(strsplit(tail(formula_split, 1), "[+]")[[1]])
    predictors1 <- predictors[predictors %in% colnames(dat)]
    #### Remove predictors that don't vary at all ####
    invalid_predictors <- predictors1[
        sapply(predictors1, function(x){length(unique(dat[[x]]))<2})
    ]
    if(length(invalid_predictors)>0){
        messager("Removing predictors that don't vary in the data:",
                 paste("\n -",invalid_predictors,collapse = ""),
                 v=verbose)
        predictors1 <- predictors1[!predictors1 %in% invalid_predictors]
    }
    #### Create formula ####
    formula1 <- formula(paste(
        formula_split[2],
        formula_split[1],
        paste(predictors1, collapse = " + ")
    ))
    return(formula1)
}
