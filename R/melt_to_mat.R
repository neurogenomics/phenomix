melt_to_mat <- function(dat,
                        formula = "trait1 ~ trait2",
                        fun.aggregate = mean, 
                        fill = 0,
                        value.var = "similarity"){
    var1 <- trimws(strsplit(formula,"~")[[1]][1])
    cmat <- data.table::dcast.data.table(data = dat,
                                         formula = formula, 
                                         fun.aggregate = fun.aggregate, 
                                         fill = fill,
                                         value.var = value.var) %>%
        tibble::column_to_rownames(var1) %>% 
        as.matrix()
    return(cmat)
}