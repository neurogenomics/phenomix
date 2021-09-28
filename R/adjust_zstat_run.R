adjust_zstat_run <- function(dat,
                             model = NULL,
                             formula,
                             verbose = TRUE,
                             ...) {
    ADJ_ZSTAT <- NULL;
    #### Select model ####
    if (is.null(model)) {
        messager("Defaulting to model stats::lm.",
            v = verbose
        )
        model <- stats::lm
    }
    #### Run ####
    messager("Adjusting dat ZSTAT.", v = verbose)
    mod <- model(
        formula = formula,
        data = dat,
        ...
    )
    #### Create adjusted predictors rowsums ####
    predictors <- names(mod$coefficients)[-1]
    adjusted_predictors <- lapply(seq(1, length(predictors)), function(i) {
        dat[[predictors[i]]] * mod$coefficients[i + 1]
    }) %>% `names<-`(predictors)
    #### Compute Adjust Z-stat ####
    dat[, ADJ_ZSTAT := dat$ZSTAT - Reduce(`+`, adjusted_predictors)]
}
